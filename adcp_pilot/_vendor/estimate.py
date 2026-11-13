"""How ADCP dockings share an allocation: how many run at once, and with
how many cores each.

By Aadhil Haq
"""

# ADCP runs each Monte Carlo replica as its own process, as many at a time as
# it has cores, so one docking stops getting faster once a single wave holds
# every replica. Measured with 20 replicas against one target:
#
#     cores   seconds
#         1       127
#         8        21
#        16        17
#        32        15
#
# A wave takes as long as its slowest replica, and a wider wave has slower
# stragglers, so beyond about eight cores extra cores barely help a single
# docking. Several dockings side by side use a large allocation better: four
# dockings took 72 s one after another on 16 cores, and 42 s as four at once
# on 4 cores each.
CORES_PER_DOCKING = 8


def plan_concurrency(cpus, n_jobs, per_docking=CORES_PER_DOCKING,
                     replicas=None):
    """(how many dockings at once, cores each).

    ``replicas`` is accepted for compatibility and ignored.
    """
    cpus = max(1, int(cpus or 1))
    n_jobs = max(1, int(n_jobs))
    workers = max(1, min(n_jobs, cpus // max(1, per_docking)))
    return workers, max(1, cpus // workers)
