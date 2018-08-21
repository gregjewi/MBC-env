import swmmAPI_v2 as sm
from pyswmm import Simulation, Links, Nodes
import numpy as np

# HoF = dict()
# with open('GA_1_outputs.csv','r') as f:
#     next(f)
#     for l in f:
#         vals = f.readline()        
#         vals = vals.rstrip().split(',')
#         key = (vals[0],vals[1])
#         HoF[key] = vals[2:] # first two values are the fitness (flood_count,total CSO)

import pickle
with open('GA_opt_20150530_GEN20.p','rb') as f:
    HoF = pickle.load(f)
        
# Run Through All the simulations:
for h in HoF:
    uparams = HoF[h][:6]
    dparams = HoF[h][6:9]
    setpts = HoF[h][9:12]
    
    swmmINP = 'NoWestside_withISDs_RainEvents_Simplifying_copy.inp'
    # swmmINP = 'NoWestside_withISDs_RainEvents_NoControl1.inp'
    # outF = 'HoF_20150530_'+h[0]+'_'+h[1]+'.out'
    outF = 'HoF_20150530_' + str(round(h[0])) + '_' + str(round(h[1])) + '.out'
    
    
    model = sm.swmmINP(swmmINP)
    model.set_dicts()
    model.set_geo_dicts()
    # model.convert(7.0)
    
    offset = -479.755 # To Detroit datum
    control = True
    
    ControlPoints = sm.make_control_points('input_files/ControlPoints.csv')
    DownstreamPoints = sm.make_downstream_points('input_files/Downstream_Init_FV.csv')

    for c,i in zip(ControlPoints,uparams):
        c.get_model_info(model)
        c.u_param = float(i)
    for d,i,j in zip(DownstreamPoints,dparams,setpts):
        d.get_model_info(model)
        d.epsilon = float(i)
        d.set_point = float(j)
        
    for d in DownstreamPoints:
        if d.d_type == 'link':
            pass
        else:
            d.dmi['q_full'] = d.max_depth
        
    with Simulation(swmmINP,outputfile = outF) as sim:
        run = sm.system(sim,offset = offset, control = control) # Offset = 0, by default
        run.timestep = model.options['ROUTING_STEP']
        nodes = Nodes(sim)
        links = Links(sim)

        run.groups = max(d.group for d in DownstreamPoints)

        for c in ControlPoints:
            c.set_vars(nodes,links)
        for d in DownstreamPoints:
            d.set_vars(nodes,links)


        # DO ONCE
        for d in DownstreamPoints:
            run.group_lookup[d.group] = d.d_type
        run.group_lookup

        n = len(ControlPoints)
        link_mask = np.zeros((1,n))
        storage_mask = np.zeros((1,n))
        for c,i in zip(ControlPoints,range(0,n)):
            if run.group_lookup[c.group] == 'link':
                link_mask[0,i] = 1.0
            elif run.group_lookup[c.group] == 'storage':
                storage_mask[0,i] = 1.0
            else:
                print('Check out', c.c_name, 'neither link nor storage')

        # make group matrix. dimensions of matrix are [ # of groups , # of Control Elements ]
        # Put 1 in row if Control Element is part of group
        groupM = np.zeros((run.groups,n))
        for c,i in zip(ControlPoints,range(0,n)):
            groupM[c.group-1,i] = 1

        # Make arrays for price calculations
        uparam = np.array([c.u_param for c in ControlPoints])
        dparam = np.array([d.epsilon for d in DownstreamPoints]) # Do Once
        setpts = np.array([d.set_point for d in DownstreamPoints]) # Do Once
        n_tanks_1 = np.sum(groupM,axis=1)+1 # Do once
        q_goals = [c.q_goal for c in ControlPoints]
        q_full = [d.dmi['q_full'] for d in DownstreamPoints]
        max_depth = np.array([d.max_depth for d in DownstreamPoints])
        set_and_full = setpts * q_full # Do once.


        print('Running Simulation... ',str(h))
        for step in sim:
            if run.control:
                ustream = np.array([c.u_var.depth / c.max_depth for c in ControlPoints])
                Uprods = uparam*ustream
                # Sum of Uparam*Ustream for each group
                np.mat(Uprods)*np.mat(groupM).transpose()

                dstream = np.array([d.d_var.depth / d.max_depth for i in DownstreamPoints])
                d_calcs = ( dstream -setpts ) * dparam # If dstream > set, P increases

                # Pareto Price for each market
                P = (np.mat(Uprods)*np.mat(groupM).transpose() + d_calcs)/ n_tanks_1
                # Give pareto price of respective group to all elements
                Pmat = np.mat(P) * np.mat(groupM)
                # Calculate individual element demands
                Pdemand = -Pmat + Uprods
                Pdemand[Pdemand < 0] = 0
                # Supply for each group
                PS = Pdemand * np.mat(groupM).transpose()

                # q_goal if all were links
                check_zero = np.divide(set_and_full, PS, out=np.zeros_like(PS), where=PS!=0)
                q_links = np.array(check_zero * groupM )* np.array( Pdemand )

                # q_goal if all were storages
                q_storage = np.array(Pdemand) * np.array( np.mat(max_depth/run.timestep) * groupM )

                # Apply masks and add together for final q_goals
                q_links = q_links * link_mask
                q_storage = q_storage * storage_mask
                q_goal = q_links + q_storage

                # Assign and send
                for q,c in zip(np.ndenumerate(q_goal),ControlPoints):
                    c.q_goal = q[1] # q is a tuple with [0] index of array and [1] the value
                    c.get_target_setting(run,nodes,links)

    run.flood_count = sum(c.flood_count for c in ControlPoints)