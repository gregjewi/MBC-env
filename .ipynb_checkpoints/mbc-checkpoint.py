import numpy as np
from pyswmm import Simulation, Links, Nodes
import pandas as pd
import pickle
import csv
import matplotlib.pyplot as plt

def new_fnctn(upstreamDict,downstreamDict,controlDict,dsKeys,price,PDemand,nodes):
    # These currently don't change per timestep. That is a possibility though
    setpts = np.array([downstreamDict[dsKeys[0]]['set_point'],downstreamDict[dsKeys[0]]['set_derivative']])
    dparam = [downstreamDict[dsKeys[0]]['epsilon'],downstreamDict[dsKeys[0]]['gamma']]
    uparam = np.array([upstreamDict[i]['uparam'] for i in upstreamDict]) 
    n_tanks = len(upstreamDict)
    
    ustream = np.array([upstreamDict[i]['ts'][-1] for i in upstreamDict])
    try:
        dstream = np.array([downstreamDict[dsKeys[0]]['ts_flow'][-1],downstreamDict[dsKeys[0]]['ts_flow'][-1]-downstreamDict[dsKeys[0]]['ts_flow'][-2]])
    except:
        dstream = np.array([downstreamDict[dsKeys[0]]['ts_flow'][-1],0])
    
    # Set pareto optimal price
    p = (sum(uparam*ustream) + sum(dparam*(dstream-setpts)))/(1 + n_tanks)

    # Initialize demand array holder
    PD = np.zeros(n_tanks)
    
    # Calculate individual Demand Prices for each upstream storage/conduit
    for i in range(0,n_tanks):
        PD[i] = max(-p + uparam[i]*ustream[i],0)
    PS = sum(PD)
    
    # Store Demand Prices for timestep in their dictionaries
    for i,j in zip(upstreamDict,range(0,len(PD))):
        upstreamDict[i]['PD'].append(PD[j])
    
    # Calculate Qi or q_goal, send to target setting
    for i,j in zip(upstreamDict,controlDict):
        if PS == 0:
            controlDict[j]['q_goal'] = 0
        else:
            # have not implemented cascasding summation.
            controlDict[j]['q_goal'] = upstreamDict[i]['PD'][-1]/PS*setpts[0]*downstreamDict[dsKeys[0]]['max_flow'] # setpts[0] assumed to be downstream flow setpoint
        
        if upstreamDict[i]['ts'][-1] == 0:
            controlDict[j]['action'] = 0.0 #Keep orifice closed
        else:
            # Send dictionary of control points. "Output" is changing the 'action', meaning target_setting.
            get_target_setting(controlDict,nodes)
    
    price.append(p)
    PDemand.append(PD)

def get_target_setting(controlDict, nodes):
    # do something to calculate target setting for each orifice/pump
    # Assign target_setting as 'action' in each control point dict
    
    # Assumption is that orifice flow only, never weir flow.
    # Need the following variables:
    #   - Depth of water in upstream node (h1)
    #   - Depth of water in downstream node (h2)
    #   - Current weir setting (before changing the setting)
    #   - Orifice Geometries
    
    g = 32.2 # gravity
    
    for i in controlDict:
        control_connects = controlDict[i]['pyswmmVar'].connections
        upstream = nodes[control_connects[0]]
        downstream = nodes[control_connects[1]]

        h1 = upstream.depth + upstream.invert_elevation
        h2 = downstream.depth + downstream.invert_elevation
        current_setting = controlDict[i]['pyswmmVar'].current_setting

        if controlDict[i]['type'] == 'SIDE': # means side orifice
            if controlDict[i]['shape'] == 'RECT_CLOSED':
                current_height = current_setting * controlDict[i]['geom1']
                h_midpt = (current_height / 2) + (upstream.invert_elevation + downstream.invert_elevation) / 2

                if h2 < h_midpt:
                    H = h1 - h_midpt
                    y = controlDict[i]['q_goal'] / (controlDict[i]['Cd'] * np.sqrt(2*g) * controlDict[i]['geom1'] * controlDict[i]['geom2'])

                    answer = np.roots([controlDict[i]['geom1']/2,h1 - upstream.invert_elevation,0,-y**2])
                    for number in answer:
                        if np.real(number) < 1.0 and np.real(number) > 0.0:
                            controlDict[i]['action'] = np.real(number)

                    # This took too long... switched to root finder above from numpy
                    # x = Symbol('x')
                    # answer = solve(controlDict[i]['geom1']/2 * x**3 + h1 * x**2 - y**2)
                    # for number in answer:
                    #     if sympy.functions.re(number) > 0.0 and sympy.functions.re(number) < 1.0:
                    #         controlDict[i]['action'] = sympy.functions.re(number)

                else:
                    H = h1 - h2
                    # print(H)
                    # if H is close to zero, don't change the target setting.
                    if np.allclose(0.0,H,atol=1e-04) or H < 0.0:
                        controlDict[i]['action'] = current_setting
                    else:
                        target_setting = controlDict[i]['q_goal'] / (controlDict[i]['Cd'] * np.sqrt(2 * g * H) * controlDict[i]['geom1'] * controlDict[i]['geom2'])
                        controlDict[i]['action'] = target_setting
            
            elif controlDict[i]['shape'] == 'CIRCULAR':
                # do something here to calculate target setting for circular orifice
                pass

        elif controlDict[i]['type'] == 'pump':
            if controlDict[i]['curve_info']['type'] == 'PUMP3':
                head = h2 - h1 # pump pushes water from low head (h1) to higher head (h2)

                # calculate q_full at given head.
                q_full = np.interp(
                    head,
                    np.array(controlDict[i]['curve_info']['x_val'],dtype='float64'),
                    np.array(controlDict[i]['curve_info']['y_val'],dtype='float64'),
                    )

                if controlDict[i]['q_goal'] == 0.0:
                    controlDict[i]['action'] = 0.0
                else:
                    controlDict[i]['action'] = q_full / controlDict[i]['q_goal']

        else:
            # do something here with the other orifice and pumps
            pass
        
        # if target setting greater than 1, only open to 1.
        controlDict[i]['action'] = min(controlDict[i]['action'], 1.0)
        # controlDict[i]['action'] = np.random.random() # See if it actually works...

def fnctn(ustream, dstream, setpts, uparam, dparam, n_tanks, action,max_flow,controlDict):
    p = (sum(uparam*ustream) + sum(dparam*(dstream-setpts)))/(1 + n_tanks)
    PD = np.zeros(n_tanks)
    for i in range(0,n_tanks):
        PD[i] = max(-p + uparam[i]*ustream[i],0)
    PS = sum(PD)
    
    for i in range(0,n_tanks):
        if PS == 0:
            Qi = 0
        else:
            # Qi = PD[i]/PS*setpts[0] # setpts[0] assumed to be downstream flow setpoint
            Qi = sum(PD[0:i+1])/PS*setpts[0]*max_flow # setpts[0] assumed to be downstream flow setpoint
        if ustream[i] == 0:
            action[i] = 0.5
        else:
            h2i = Qi/(1.0*1*np.sqrt(2*9.81*ustream[i]))
#             h2i = Qi/(0.61*1*np.sqrt(2*9.81*ustream[i]))
            action[i] = max(min(h2i/2,1.0),0.0)
#         if ustream[i] > 0.95:
#             action[i] = 1.0

    return p, PD, PS, action

def run_control_sim(controlDict,upstreamDict,downstreamDict,dsKeys,swmmINP,performanceDict):
    with Simulation(swmmINP) as sim:
        freq = '12s'
        timesteps = pd.date_range(sim.start_time,sim.end_time,freq=freq)
        price = []
        PDemand = []
        
        nodes = Nodes(sim)
        links = Links(sim)
        
        for point in controlDict:
            controlDict[point]['pyswmmVar'] = links[point]
        for point in upstreamDict:
            try:
                upstreamDict[point]['pyswmmVar'] = links[point]
            except:
                upstreamDict[point]['pyswmmVar'] = nodes[point]
        for point in downstreamDict:
            downstreamDict[point]['pyswmmVar'] = links[point]
        for point in performanceDict:
            if performanceDict[point]['type'] == 'orifice' or performanceDict[point] == 'junction' or performanceDict[point] == 'storage':
                performanceDict[point]['pyswmmVar'] = nodes[point]
            elif performanceDict[point]['type'] == 'link':
                performanceDict[point]['pyswmmVar'] = links[point]

        
        # These currently don't change per timestep. That is a possibility though
        setpts = np.array([downstreamDict[dsKeys[0]]['set_point'],downstreamDict[dsKeys[0]]['set_derivative']])
        dparam = [downstreamDict[dsKeys[0]]['epsilon'],downstreamDict[dsKeys[0]]['gamma']]
        uparam = np.array([upstreamDict[i]['uparam'] for i in upstreamDict]) #for now make it a scalar... can diversify if needed once we want to prioritize upstream links above others
        n_tanks = len(upstreamDict)

        print('running simulation...')
        for step in sim:
            for point in upstreamDict:
                upstreamDict[point]['ts'].append(upstreamDict[point]['pyswmmVar'].depth / upstreamDict[point]['max_depth'])
            for point in downstreamDict:
                downstreamDict[point]['ts_flow'].append(downstreamDict[point]['pyswmmVar'].flow / downstreamDict[point]['max_flow'])
                downstreamDict[point]['ts_depth'].append(downstreamDict[point]['pyswmmVar'].depth / downstreamDict[point]['max_depth'])
            for point in performanceDict:
                try:
                    performanceDict[point]['ts_flow'].append(performanceDict[point]['pyswmmVar'].flow)
                except:
                    pass

                try:
                    performanceDict[point]['ts_flow'].append(performanceDict[point]['pyswmmVar'].total_inflow)
                except:
                    pass
            
            ustream = np.array([upstreamDict[i]['ts'][-1] for i in upstreamDict])
            try:
                dstream = np.array([downstreamDict[dsKeys[0]]['ts_flow'][-1],downstreamDict[dsKeys[0]]['ts_flow'][-1]-downstreamDict[dsKeys[0]]['ts_flow'][-2]])
            except:
                dstream = np.array([downstreamDict[dsKeys[0]]['ts_flow'][-1],0])
            
            
            action = [controlDict[i]['action'] for i in controlDict]
            p, PD, PS, action = fnctn(ustream, dstream, setpts, uparam, dparam, n_tanks, action,downstreamDict[dsKeys[0]]['max_flow'],controlDict)
            
            # try new function for mbc.
            new_fnctn(upstreamDict,downstreamDict,controlDict, dsKeys, price, PDemand, nodes)
            # print([controlDict[i]['action'] for i in controlDict])
            
            # for point,i in zip(controlDict,range(0,len(action))):
            for point in controlDict:
                controlDict[point]['pyswmmVar'].target_setting = controlDict[point]['action']
                # controlDict[point]['actionTS'].append(action[i])
                controlDict[point]['actionTS'].append(controlDict[point]['action'])
                
                
    print('done ...')
    return timesteps,price,PDemand

def run_no_control_sim(controlDict,upstreamDict,downstreamDict,dsKeys,swmmINP,performanceDict):
    with Simulation(swmmINP) as sim:
        freq = '12s'
        timesteps = pd.date_range(sim.start_time,sim.end_time,freq=freq)
        price = []
        PDemand = []
        
        nodes = Nodes(sim)
        links = Links(sim)
        
        for point in controlDict:
            controlDict[point]['pyswmmVar'] = links[point]
        for point in upstreamDict:
            upstreamDict[point]['pyswmmVar'] = links[point]
        for point in downstreamDict:
            downstreamDict[point]['pyswmmVar'] = links[point]
        for point in performanceDict:
            if performanceDict[point]['type'] == 'orifice' or performanceDict[point] == 'junction' or performanceDict[point] == 'storage':
                performanceDict[point]['pyswmmVar'] = nodes[point]
            elif performanceDict[point]['type'] == 'link':
                performanceDict[point]['pyswmmVar'] = links[point]

        print('running simulation...')
        for step in sim:
            for point in controlDict:
                controlDict[point]['pyswmmVar'].target_setting = 1.0
            for point in upstreamDict:
                upstreamDict[point]['ts'].append(upstreamDict[point]['pyswmmVar'].depth / upstreamDict[point]['max_depth'])
            for point in downstreamDict:
                downstreamDict[point]['ts_flow'].append(downstreamDict[point]['pyswmmVar'].flow / downstreamDict[point]['max_flow'])
                downstreamDict[point]['ts_depth'].append(downstreamDict[point]['pyswmmVar'].depth / downstreamDict[point]['max_depth'])
            for point in performanceDict:
                try:
                    performanceDict[point]['ts_flow'].append(performanceDict[point]['pyswmmVar'].flow)
                except:
                    pass

                try:
                    performanceDict[point]['ts_flow'].append(performanceDict[point]['pyswmmVar'].total_inflow)
                except:
                    pass
    print('done ...')

def save_that_shit(saveDict):
    # Save outputs as pickle. Pickle can't save the pyswmm objects so we need to get rid of those.
    for key in saveDict:
        if key == 'controlDict' or key == 'downstreamDict' or key == 'upstreamDict' or key=='performanceDict':
            for var in saveDict[key]:
                try:
                    saveDict[key][var].pop('pyswmmVar')
                except:
                    pass

    with open(saveDict['pickleOut'],'wb') as f:
        pickle.dump(saveDict,f)

    metadata = []
    metadata.append(saveDict['pickleOut'])
    metadata.append(saveDict['file'])
    metadata.append(saveDict['control'])
    varL = list(saveDict.keys())

    if 'upstreamDict' in varL:
        for i in saveDict['upstreamDict']:
            metadata = metadata + [saveDict['upstreamDict'][i]['DScp'],saveDict['upstreamDict'][i]['beta'],saveDict['upstreamDict'][i]['uparam']]
    else:
        for i in saveDict['upstreamDict']:
            metadata = metadata + ['na']

    if 'downstreamDict' in varL:
        for i in saveDict['downstreamDict']:
            metadata = metadata + [i,saveDict['downstreamDict'][i]['epsilon'],saveDict['downstreamDict'][i]['gamma'],saveDict['downstreamDict'][i]['max_flow'],saveDict['downstreamDict'][i]['max_depth'],saveDict['downstreamDict'][i]['set_point'],saveDict['downstreamDict'][i]['set_derivative']]
    else:
        for i in saveDict['downstreamDict']:
            metadata = metadata + ['na']

    if 'notes' in varL:
        for i in saveDict['notes']:
            metadata = metadata + [i]
        else:
            pass

    with open(saveDict['metaCSV'],'a+',newline='') as f:
        writer = csv.writer(f)
        writer.writerow(metadata)

def viz_shit(f_nocontrol,f_control,figname):
    # Unpack the pickles
    with open(f_control,'rb') as fC, open(f_nocontrol,'rb') as fNC:
        data_control = pickle.load(fC)
        data_nocontrol = pickle.load(fNC)

    fig,axarr = plt.subplots(2,2, figsize=(12,8))

    # Plot Control Points onto the figures
    for p1,p2 in zip(data_control['upstreamDict'],data_nocontrol['upstreamDict']):
        axarr[0,0].plot(data_control['upstreamDict'][p1]['ts'])
        axarr[0,0].plot(data_nocontrol['upstreamDict'][p2]['ts'],linestyle=':')
    axarr[0,0].set_title('Upstream Depth')
    
    # make legend
    leg = []
    leg = leg + [data_control['upstreamDict'][i]['DScp'] for i in data_control['upstreamDict']]
    leg = leg + [data_nocontrol['upstreamDict'][i]['DScp'] for i in data_nocontrol['upstreamDict']]
    axarr[0,0].legend(leg,title = 'Control__, No Control ..', ncol=2,loc=4)
    axarr[0,0].set_ylabel('Normalized Depth')
    # leg_control = plt.legend([data_control['upstreamDict'][i]['DScp'] for i in data_control['upstreamDict']])
    # axarr[0,0].add_artist(leg_control)
    # leg_nocontrol = [data_nocontrol['upstreamDict'][i]['DScp'] for i in data_nocontrol['upstreamDict']]
    # axarr[0,0].legend(leg_nocontrol)
    # plt.show()

    for p1,p2 in zip(data_control['downstreamDict'],data_nocontrol['downstreamDict']):
        axarr[0,1].plot(data_control['downstreamDict'][p1]['ts_flow'])
        axarr[0,1].plot(data_nocontrol['downstreamDict'][p2]['ts_flow'],color='k',linestyle=':')
    axarr[0,1].set_title('Downstream Flow')
    axarr[0,1].legend(['Control', 'No Control'])
    axarr[0,1].set_ylabel('Normalized Flow')
    # plt.show()

    for p1,p2 in zip(data_control['downstreamDict'],data_nocontrol['downstreamDict']):
        axarr[1,0].plot(data_control['downstreamDict'][p1]['ts_depth'])
        axarr[1,0].plot(data_nocontrol['downstreamDict'][p2]['ts_depth'],color='k',linestyle=':')
    axarr[1,0].set_title('Downstream Depth')
    axarr[1,0].legend(['Control','No Control'])
    axarr[1,0].set_ylabel('Noramlized Depth')

    # axarr[1,1].plot(data_control['performanceDict']['1002'])
    # axarr[1,1].plot(data_nocontrol['performanceDict']['1002'],color='k',linestyle=':')
    plt.show()

    # plt.savefig(figname)
    plt.close()

    # plt.plot()
    # plt.plot(data_control['performanceDict']['1002'])
    # plt.plot(data_nocontrol['performanceDict']['1002'],color='k',linestyle=':')
    # # axarr.set_title('Inflow to WWTP from DRI')
    # # axarr.legend(['Control','No Control'])
    # # axarr.set_ylabel('Flow [cfs]')
    # plt.show()
    # plt.close()

    # return fig