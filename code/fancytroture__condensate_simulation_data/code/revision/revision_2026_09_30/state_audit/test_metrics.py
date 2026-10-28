"""Numerical invariants for missing-state and mask sensitivity calculations."""
import numpy as np
from audit_states import metric, distance


def main():
    t = np.array([[.7,.2,.1],[.1,.8,.1],[.2,.2,.6]])
    p = np.array([.4,.4,.2])
    common = np.array([True,True,False])
    one = np.ones((3,3))
    for definition in ("unfiltered_display_union","common_submatrix","common_conditional","common_plus_exit"):
        assert metric(definition,t,t,p,p,common,one,one) == 0
    changed_mask = one.copy()
    changed_mask[0,1] = 0
    assert metric("submitted_edge_mask",t,t,p,p,common,one,changed_mask) > 0
    u = t.copy()
    u[0] = [.6,.2,.2]
    a = np.array([[.7,.2,.1],[.1,.8,.1]])
    b = np.array([[.6,.2,.2],[.1,.8,.1]])
    expected = distance(a,b,np.array([.5,.5]),np.array([.5,.5]))
    actual = metric("common_plus_exit",t,u,p,p,common,one,one)
    assert np.isclose(expected,actual)
    # Restricting to C and conditioning on staying in C answers a different question.
    assert not np.isclose(metric("common_conditional",t,u,p,p,common,one,one),actual)
    perm = np.array([2,0,1])
    assert np.isclose(actual,metric("common_plus_exit",t[np.ix_(perm,perm)],u[np.ix_(perm,perm)],
                                   p[perm],p[perm],common[perm],one,one))
    batched = metric("common_plus_exit",np.stack([t,t]),np.stack([u,u]),
                    np.stack([p,p]),np.stack([p,p]),common,one,one)
    assert np.allclose(batched,actual)
    print("Passed: identity, edge-mask artifact, escape mass, conditioning distinction, permutation, batch.")


if __name__ == "__main__":
    main()
