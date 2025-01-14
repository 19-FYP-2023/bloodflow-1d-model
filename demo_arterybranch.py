import sys

import numpy as np
import configparser

import arteryfe as af


def main(config_location):
    """Read config-file.
    Run the necessary functions to compute the solution.
    :param string config_location: Location of config file
    """
    param = af.ParamParser(config_location)

    # Constructor parameters
    no_of_arteries = param.param['no_of_arteries']
    rc = param.param['rc']
    qc = param.param['qc']
    k1 = param.param['k1']
    k2 = param.param['k2']
    k3 = param.param['k3']
    rho = param.param['rho']
    nu = param.param['nu']
    p0 = param.param['p0']
    geometric_data_location = param.param['geometric_data_location']

    # Geometry parameters
    Nt = param.geo['Nt']
    Nx = param.geo['Nx']
    N_cycles = param.geo['N_cycles']

    # Solution parameters
    inlet_flow_location = param.solution['inlet_flow_location']
    output_location = param.solution['output_location']
    theta = param.solution['theta']
    Nt_store = param.solution['Nt_store']
    N_cycles_store = param.solution['N_cycles_store']
    store_area = param.solution['store_area']
    store_pressure = param.solution['store_pressure']

    # Import inlet flow data
    T, q_ins = af.read_inlet(inlet_flow_location, Nt)
    
    # nondimensionalize inlet flow
    q_ins = q_ins/ qc
    print("q_ins: ", q_ins[0])
    T = T*qc/rc**3

    # Create artery network
    an = af.ArteryNetwork(no_of_arteries, rc, qc, k1, k2, k3, rho, nu, p0, geometric_data_location)
    
    an.define_geometry(Nx, Nt, T, N_cycles)
    an.define_solution(output_location, q_ins[0], theta)

    # Solve problem and store data
    an.solve(q_ins, Nt_store, N_cycles_store, store_area, store_pressure)


if __name__ == '__main__':
    main(sys.argv[1])
