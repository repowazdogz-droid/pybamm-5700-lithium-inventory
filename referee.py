"""Independent finite-volume inventory and analytic fault controls."""
import json
import numpy as np


def inventory(area, electrolyte, particles):
    total = None
    for concentration, widths, fraction in electrolyte:
        value = area * fraction * np.einsum('i,it->t', widths, concentration)
        total = value if total is None else total + value
    for concentration, widths, edges, fraction in particles:
        weights = np.diff(edges**3) / (edges[-1]**3 - edges[0]**3)
        total += area * fraction * np.einsum('i,j,jit->t', widths, weights, concentration)
    return total


def inventory_loop(area, electrolyte, particles):
    nt = electrolyte[0][0].shape[-1]
    answer = np.zeros(nt)
    for c, dx, eps in electrolyte:
        for i in range(len(dx)):
            answer += area * eps * dx[i] * c[i, :]
    for c, dx, r, eps in particles:
        volume = 4 * np.pi * (r[-1]**3-r[0]**3) / 3
        for i in range(len(dx)):
            for j in range(len(r)-1):
                shell = 4 * np.pi * (r[j+1]**3-r[j]**3) / 3
                answer += area * eps * dx[i] * shell / volume * c[j,i,:]
    return answer


def residual(n, exchange, source=0):
    return n-n[0]-exchange-source


def controls():
    import pint
    import sympy as sp
    u=pint.UnitRegistry()
    assert (u.meter**2*u.meter*u.mole/u.meter**3).dimensionality == u.mole.dimensionality
    assert (u.ampere*u.second/(u.coulomb/u.mole)).dimensionality == u.mole.dimensionality
    assert (u.ampere/u.meter**3/(u.coulomb/u.mole)).dimensionality == (u.mole/u.meter**3/u.second).dimensionality
    assert (u.mole/(u.coulomb/u.mole)).dimensionality != u.mole.dimensionality
    x=sp.Symbol('x'); t=sp.Function('t')(x); i=sp.Function('i')(x)
    assert sp.simplify(sp.diff(t*i,x)-t*sp.diff(i,x)-i*sp.diff(t,x))==0
    times=np.linspace(0,1,101); dx=np.array([0.2,0.8]); r=np.array([0.,0.2,0.7,1.])
    area=2.;eps_e=0.4;eps_s=0.6
    ce=np.tile(3+times,(2,1));cs=np.tile(8-eps_e/eps_s*times,(3,2,1))
    e=[(ce,dx,eps_e)];p=[(cs,dx,r,eps_s)]
    n=inventory(area,e,p);n_loop=inventory_loop(area,e,p)
    np.testing.assert_allclose(n,n_loop,rtol=1e-14,atol=1e-14)
    q=0.3*times; s=0.2*times**2; open_n=n+q+s
    values={
      'closed_transfer':residual(n,0),
      'boundary_and_source_conserving':residual(open_n,q,s),
      'altered_inventory':residual(open_n+0.02*times,q,s),
      'wrong_boundary_sign':residual(open_n,-q,s),
      'omitted_source':residual(open_n,q),
      'wrong_volume_factor':residual(open_n*1.01,q,s),
      'synthetic_imbalance':residual(open_n+0.01*np.sin(np.pi*times),q,s)
    }
    maxima={k:float(np.max(np.abs(v))) for k,v in values.items()}
    for k,v in maxima.items():
        assert (v<1e-12) if k in ['closed_transfer','boundary_and_source_conserving'] else (v>1e-4), (k,v)
    return {'passed':True,'max_residual_mol':maxima,'loop_disagreement':float(max(abs(n-n_loop))), 'units_and_product_rule':'passed'}


if __name__=='__main__':
    print(json.dumps(controls(),indent=2))
