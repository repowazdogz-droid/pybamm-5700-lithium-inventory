"""Explicit PyBaMM API contract and independent lithium measurements."""
import json
import pathlib
import numpy as np
import pybamm
from scipy.integrate import cumulative_trapezoid, simpson
from referee import inventory, inventory_loop

ROOT = pathlib.Path(__file__).resolve().parent
FIELDS = {
    'negative electrode': 'Negative electrolyte concentration [mol.m-3]',
    'separator': 'Separator electrolyte concentration [mol.m-3]',
    'positive electrode': 'Positive electrolyte concentration [mol.m-3]',
}
GEOMETRY = {
    'negative electrode': ('Negative electrode porosity', 'Negative electrode thickness [m]'),
    'separator': ('Separator porosity', 'Separator thickness [m]'),
    'positive electrode': ('Positive electrode porosity', 'Positive electrode thickness [m]'),
}
PARTICLES = {
    'negative': 'Negative particle concentration [mol.m-3]',
    'positive': 'Positive particle concentration [mol.m-3]',
}
STEPS = [
    'Rest for 10 minutes (5 seconds period)',
    'Discharge at 1C for 20 minutes (5 seconds period)',
    'Rest for 10 minutes (5 seconds period)',
    'Charge at 1C for 20 minutes (5 seconds period)',
    'Rest for 20 minutes (5 seconds period)',
]


class MeasurementError(RuntimeError):
    pass


def require(condition, message):
    if not condition:
        raise MeasurementError(message)


def prepare(case, nx, nr, tol, smoke=False):
    basic = case.startswith('basic')
    half = case == 'half'
    model = pybamm.lithium_ion.BasicDFN() if basic else pybamm.lithium_ion.DFN({'working electrode': 'positive'} if half else {})
    parameters = pybamm.ParameterValues('Xu2019' if half else 'ORegan2022')
    if not half:
        parameters['Negative particle diffusivity [m2.s-1]'] = 3.3e-14
    if case == 'basic_constant':
        parameters['Cation transference number'] = parameters['Cation transference number'](parameters['Initial concentration in electrolyte [mol.m-3]'], parameters['Ambient temperature [K]'])
    domains = ['separator', 'positive electrode'] if half else list(FIELDS)
    electrodes = ['positive'] if half else list(PARTICLES)
    used = [FIELDS[d] for d in domains] + [PARTICLES[e] for e in electrodes] + ['Current [A]', 'Discharge capacity [A.h]']
    for name in used:
        require(name in model.variables, 'missing API field: ' + name)
    rhs_by_name = {k.name: value for k, value in model.rhs.items()}
    for electrode in electrodes:
        name = PARTICLES[electrode]
        require(name in rhs_by_name, 'missing particle RHS: ' + name)
        model.variables['Study rate ' + name] = rhs_by_name[name]
    rhs_name = 'Electrolyte concentration [mol.m-3]' if basic else 'Porosity times concentration [mol.m-3]'
    require(rhs_name in rhs_by_name, 'missing electrolyte RHS')
    model.variables['Study electrolyte rate'] = rhs_by_name[rhs_name]
    if basic:
        supplied = pybamm.Scalar(0)
        for domain in domains:
            eps, length = GEOMETRY[domain]
            supplied += pybamm.x_average(model.variables[FIELDS[domain]] * pybamm.Parameter(eps)) * pybamm.Parameter(length)
        for electrode in electrodes:
            title = electrode.capitalize()
            c = model.variables[PARTICLES[electrode]]
            supplied += pybamm.x_average(pybamm.Parameter(title + ' electrode active material volume fraction') * pybamm.r_average(c)) * pybamm.Parameter(title + ' electrode thickness [m]')
            phase = getattr(model.param, electrode[0]).prim
            j_over_f = -model.boundary_conditions[c]['right'][0] * pybamm.surf(phase.D(c, model.param.T_init))
            aj_over_f = 3 * phase.epsilon_s_av / phase.R_typ * j_over_f
            ce = model.variables[FIELDS[electrode + ' electrode']]
            model.variables['Study predictor ' + electrode] = -model.param.t_plus(ce, model.param.T_init) * aj_over_f
        model.variables['Study supplied inventory'] = supplied / pybamm.constants.F.value * model.param.A_cc
    else:
        model.variables['Study supplied inventory'] = model.variables['Total lithium [mol]']
    solver = pybamm.IDAKLUSolver(rtol=tol, atol=tol, root_method='casadi', root_tol=1e-10, options={'num_threads': 1})
    experiment = pybamm.Experiment(['Rest for 5 seconds (1 second period)', 'Discharge at C/20 for 5 seconds (1 second period)'] if smoke else STEPS)
    simulation = pybamm.Simulation(model, parameter_values=parameters, experiment=experiment, solver=solver, var_pts={'x_n': nx, 'x_s': nx, 'x_p': nx, 'r_n': nr, 'r_p': nr})
    return simulation, parameters, domains, electrodes, used


def extract(simulation, solution, parameters, domains, electrodes, case, destination, smoke=False):
    """Same extraction path in preflight and scientific runs; no conservation threshold here."""
    area = float(parameters['Electrode width [m]'] * parameters['Electrode height [m]'] * parameters['Number of electrodes connected in parallel to make a cell'])
    require(np.isfinite(area) and area > 0, 'invalid effective area')
    faraday = float(pybamm.constants.F.value)
    dx = {d: np.diff(simulation.mesh[d].edges) for d in domains}
    radial = {e: simulation.mesh[e+' particle'].edges for e in electrodes}
    for domain in domains:
        require(np.all(dx[domain] > 0), 'invalid cell widths')
        np.testing.assert_allclose(sum(dx[domain]), parameters[GEOMETRY[domain][1]], rtol=1e-12)
    for edges in radial.values():
        require(edges[0] == 0 and np.all(np.diff(edges) > 0), 'invalid spherical shell geometry')
    steps = [s for cycle in solution.cycles if cycle is not None for s in cycle.steps]
    tables = []
    external_before = 0.
    predicted_before = 0.
    simpson_total = 0.
    n0 = None
    used_shapes = {}
    metadata = {'area_m2': area, 'F_C_mol': faraday, 'field_names': dict(FIELDS), 'geometry_parameters': GEOMETRY, 'step_count': len(steps)}
    destination.mkdir(parents=True, exist_ok=True)
    for index, step in enumerate(steps):
        native_t = np.asarray(step.t)
        require(np.all(np.isfinite(native_t)) and np.all(np.diff(native_t)>0), 'nonmonotone native time')
        period = 1. if smoke else 5.
        times = np.linspace(native_t[0], native_t[-1], int(round((native_t[-1]-native_t[0])/period))+1)
        require(len(times)>1 and np.all(np.diff(times)>0), 'invalid measurement time')
        nt = len(times)
        arrays = {'time_s': times, 'native_time_s': native_t, 'area_m2': np.array(area)}
        def field(name, shape):
            coordinates = {}
            for domain in domains:
                if name == FIELDS[domain]:
                    coordinates = {'x': simulation.mesh[domain].nodes}
            for electrode in electrodes:
                if name in [PARTICLES[electrode], 'Study rate '+PARTICLES[electrode]]:
                    coordinates = {'r': simulation.mesh[electrode+' particle'].nodes, 'x': simulation.mesh[electrode+' electrode'].nodes}
                elif name == 'Study predictor '+electrode:
                    coordinates = {'x': simulation.mesh[electrode+' electrode'].nodes}
            if name == 'Study electrolyte rate':
                coordinates = {'x': simulation.mesh[domains].nodes}
            elif name == 'Electrolyte flux [mol.m-2.s-1]':
                coordinates = {'x': simulation.mesh[domains].edges}
            data = np.asarray(step[name](times, **coordinates), dtype=float)
            require(data.shape == shape, f'{name}: shape {data.shape}, expected {shape}')
            require(np.all(np.isfinite(data)), 'nonfinite field: ' + name)
            used_shapes[name] = list(data.shape[:-1]) + ['time']
            return data
        electrolyte = []
        particles = []
        particle_rates = []
        for domain in domains:
            ce = field(FIELDS[domain], (len(dx[domain]), nt))
            fraction = float(parameters[GEOMETRY[domain][0]])
            require(0 < fraction < 1, 'invalid porosity')
            electrolyte.append((ce, dx[domain], fraction))
            arrays['ce_' + domain] = ce
            arrays['dx_' + domain] = dx[domain]
            arrays['eps_' + domain] = np.array(fraction)
        for electrode in electrodes:
            nx = len(dx[electrode+' electrode']); nr = len(radial[electrode])-1
            cs = field(PARTICLES[electrode], (nr,nx,nt))
            rate = field('Study rate ' + PARTICLES[electrode], (nr,nx,nt))
            fraction = float(parameters[electrode.capitalize()+' electrode active material volume fraction'])
            require(0 < fraction < 1, 'invalid active fraction')
            args = (dx[electrode+' electrode'], radial[electrode], fraction)
            particles.append((cs, *args)); particle_rates.append((rate, *args))
            arrays['cs_' + electrode] = cs
            arrays['cs_rate_' + electrode] = rate
            arrays['r_' + electrode] = radial[electrode]
            arrays['eps_s_' + electrode] = np.array(fraction)
        current = field('Current [A]', (nt,))
        capacity = field('Discharge capacity [A.h]', (nt,))
        n = inventory(area, electrolyte, particles)
        loop = inventory_loop(area, electrolyte, particles)
        require(np.all(np.isfinite(n)) and np.all(n>0), 'invalid independent inventory')
        np.testing.assert_allclose(n, loop, rtol=1e-12, atol=1e-16)
        arrays.update(inventory_mol=n, loop_inventory_mol=loop, current_A=current, capacity_Ah=capacity)
        # Checkpoint core measurements before optional diagnostic/rate processing.
        np.savez_compressed(destination/f'step{index}-core.npz', **arrays)
        if n0 is None:
            n0 = float(n[0])
        external_rate = current / faraday if case=='half' else np.zeros(nt)
        q = external_before + cumulative_trapezoid(external_rate, times, initial=0)
        external_before = float(q[-1])
        erate = field('Study electrolyte rate', (sum(len(v) for v in dx.values()), nt))
        erates = []; offset=0
        for domain, (_, widths, fraction) in zip(domains, electrolyte):
            erates.append((erate[offset:offset+len(widths)], widths, fraction if case.startswith('basic') else 1.))
            offset += len(widths)
        dn = inventory(area, erates, particle_rates)
        arrays['electrolyte_rate'] = erate
        if case.startswith('basic'):
            predicted_rate = np.zeros(nt)
            for electrode in electrodes:
                vals = field('Study predictor '+electrode, (len(dx[electrode+' electrode']), nt))
                predicted_rate += area*np.einsum('i,it->t', dx[electrode+' electrode'], vals)
                arrays['predictor_' + electrode] = vals
        else:
            flux = field('Electrolyte flux [mol.m-2.s-1]', (sum(len(v) for v in dx.values())+1, nt))
            predicted_rate = area*(flux[0]-flux[-1])-external_rate
            arrays['electrolyte_flux'] = flux
        integrated_predicted = predicted_before + cumulative_trapezoid(predicted_rate, times, initial=0)
        predicted_before = float(integrated_predicted[-1])
        simpson_total += float(simpson(predicted_rate, x=times))
        supplied = field('Study supplied inventory', (nt,))
        residual = n-n0-q
        table = np.column_stack([np.full(nt,index), times,n,q,residual,dn,external_rate,dn-external_rate,predicted_rate,integrated_predicted,loop,supplied,current,capacity])
        require(np.all(np.isfinite(table)), 'nonfinite final table')
        arrays['table'] = table
        np.savez_compressed(destination/f'step{index}.npz', **arrays)
        tables.append(table)
    table = np.concatenate(tables)
    require(np.all(np.diff(table[:,1]) >= 0), 'time went backwards between steps')
    np.savetxt(destination/'series.csv', table, delimiter=',', comments='', header='step,time_s,inventory_mol,input_mol,residual_mol,inventory_rate_mol_s,input_rate_mol_s,balance_rate_mol_s,predicted_rate_mol_s,predicted_integral_mol,loop_inventory_mol,supplied_inventory,current_A,capacity_Ah')
    metadata.update(field_shapes=used_shapes, initial_inventory_mol=n0, final_time_s=float(table[-1,1]), max_abs_residual_mol=float(max(abs(table[:,4]))), relative_max=float(max(abs(table[:,4]))/n0), final_residual_mol=float(table[-1,4]), max_loop_disagreement_mol=float(max(abs(table[:,2]-table[:,10]))), max_rate_mol_s=float(max(abs(table[:,7]))), max_rate_prediction_error_mol_s=float(max(abs(table[:,7]-table[:,8]))), final_prediction_mol=float(table[-1,9]), prediction_simpson_mol=simpson_total, quadrature_disagreement_mol=abs(simpson_total-table[-1,9]), max_supplied_inventory_disagreement_mol=float(max(abs(table[:,11]-table[:,2]))), max_capacity_input_disagreement_mol=float(max(abs(table[:,3]-table[:,13]*3600/faraday))) if case=='half' else 0.)
    (destination/'measurement.json').write_text(json.dumps(metadata,indent=2))
    return metadata
