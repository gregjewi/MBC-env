import numpy as np
from pyswmm import Simulation, Links, Nodes
import pandas as pd
import pickle
import csv
import matplotlib.pyplot as plt
import winsound

def new_fnctn(upstreamDict,downstreamDict,controlDict,dsKeys,price,PDemand,nodes,timestep,units):
    
    # These currently don't change per timestep. That is a possibility though
    setpts = np.array([downstreamDict[dsKeys[0]]['set_point'],downstreamDict[dsKeys[0]]['set_derivative']])
    dparam = [downstreamDict[dsKeys[0]]['epsilon'],downstreamDict[dsKeys[0]]['gamma']]
    uparam = np.array([upstreamDict[i]['uparam'] for i in upstreamDict]) 
    n_tanks = len(upstreamDict)
    
    ustream = np.array([upstreamDict[i]['ts'][-1] for i in upstreamDict])
    
    
    # ---- Section could be subfunction
    
    # if downstream set point is flow
    if downstreamDict[dsKeys[0]] == 'link':
        try:
            dstream = np.array([downstreamDict[dsKeys[0]]['ts_flow'][-1],downstreamDict[dsKeys[0]]['ts_flow'][-1]-downstreamDict[dsKeys[0]]['ts_flow'][-2]])
        except:
            dstream = np.array([downstreamDict[dsKeys[0]]['ts_flow'][-1],0])

    # else downstream set point in depth/storage... need to change if make DS setpoint an outfall.
    else:
        try:
            dstream = np.array([downstreamDict[dsKeys[0]]['ts_depth'][-1],downstreamDict[dsKeys[0]]['ts_depth'][-1]-downstreamDict[dsKeys[0]]['ts_depth'][-2]])
        except:
            dstream = np.array([downstreamDict[dsKeys[0]]['ts_depth'][-1],0])
    
    # ---- end section
    
    
    
    # Set pareto optimal price
    p = (sum(uparam*ustream) + sum(dparam*(setpts-dstream)))/(1 + n_tanks)

    # Initialize demand array holder
    PD = np.zeros(n_tanks)
    
    # Calculate individual Demand Prices for each upstream storage/conduit
    for i in range(0,n_tanks):
        # Negative price doesn't mean much... that's why make it 0.
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
            if downstreamDict[dsKeys[0]]['type'] == 'link':
                # have not implemented cascasding summation for ISDs in series, etc.
                controlDict[j]['q_goal'] = upstreamDict[i]['PD'][-1]/PS*setpts[0]*downstreamDict[dsKeys[0]]['max_flow'] # setpts[0] assumed to be downstream flow/depth setpoint
            elif downstreamDict[dsKeys[0]]['type'] == 'storage':
                # controlDict[j]['q_goal'] = upstreamDict[i]['PD'][-1]*downstreamDict[dsKeys[0]]['total_storage']
                controlDict[j]['q_goal'] = upstreamDict[i]['PD'][-1]*downstreamDict[dsKeys[0]]['total_storage']/timestep
            elif downstreamDict[dsKeys[0]]['type'] == 'outfall':
                # do something here
                pass
            else:
                pass
        
        if upstreamDict[i]['ts'][-1] == 0:
            controlDict[j]['action'] = 0.0 #Keep orifice closed
            
    # Send dictionary of control points. "Output" is changing the 'action', meaning target_setting.
    get_target_setting(controlDict,nodes, units)
    
    price.append(p)
    PDemand.append(PD)

def get_target_setting(controlDict, nodes, units):
    # do something to calculate target setting for each orifice/pump
    # Assign target_setting as 'action' in each control point dict
    
    # Assumption is that orifice flow only, never weir flow.
    # Need the following variables:
    #   - Depth of water in upstream node (h1)
    #   - Depth of water in downstream node (h2)
    #   - Current weir setting (before changing the setting)
    #   - Orifice Geometries

    # Pumps:
    # So far only pumps using type 3 pump curves have been assimilated into
    # the program. TBD on the others.
    
    # units should be a global variable.
    if units == 'US':
        g = 32.2 # gravity
    else:
        g = 9.81 # gravity
    

    for i in controlDict:
        control_connects = controlDict[i]['pyswmmVar'].connections
        upstream = nodes[control_connects[0]]
        downstream = nodes[control_connects[1]]

        h1 = upstream.depth + upstream.invert_elevation
        h2 = downstream.depth + downstream.invert_elevation
        
        current_setting = controlDict[i]['pyswmmVar'].current_setting # current_setting == hcrown

        pump = False
        if controlDict[i]['type'] == 'pump':
            pump = True
            

        if not pump:
            # Orifice Stuff
            
            current_height = current_setting * controlDict[i]['geom1'] # current_height == hcrown
            h_midpt = (current_height / 2) + (upstream.invert_elevation + downstream.invert_elevation) / 2
            hcrest = upstream.invert_elevation + controlDict[i]['offset']
            
            # inlet submergence
            if h1 < current_height:
                f = (h1 - hcrest) / (current_height - hcrest) # weir equation
            else:
                f = 1.0 # submerged.
    
            # which head to use
            if f < 1.0:
                H = h1 - hcrest
            elif h2 < h_midpt:
                H = h1 - h_midpt
            else:
                H = h1 - h2
            
            # USE CALCULATED HEAD AND DESIRED FLOW TO DETERMINE GATE OPENING ACTION
            
            # no head at orifice
            if H < 0.1 or f <= 0.0:
                controlDict[i]['action'] = 0.0
                # print('Head too small')
            elif h2 > h1:
                controlDict[i]['action'] = 0.0
                print('Backward flow condition, orifice closed')
            
            # Weir Flow
            elif (f < 1.0 and H > 0.1):
                A_open = controlDict[i]['q_goal'] / ( controlDict[i]['Cd'] * np.sqrt(2*g*H) * (2.0/3.0) )
                
                if controlDict[i]['shape'] == 'CIRCULAR':
                    print("Circular does not work yet. Action = 0.0")
                    controlDict[i]['action'] = 0.0
                else:
                    A_ratio = A_open / ( controlDict[i]['geom1'] * controlDict[i]['geom2'] )
                    controlDict[i]['action'] = A_ratio
            
            # True orifice flow
            else:
                # since q = Cd * A_open * sqrt( 2 g H )
                A_open = controlDict[i]['q_goal'] / ( controlDict[i]['Cd'] * np.sqrt(2*g*H) )
                
                if controlDict[i]['shape'] == 'CIRCULAR':
                    print("Circular does not work yet. Action = 0.0")
                    controlDict[i]['action'] = 0.0
                else:
                    A_ratio = A_open / ( controlDict[i]['geom1'] * controlDict[i]['geom2'] )
                    controlDict[i]['action'] = A_ratio
                    
    
        # Pump is true
        else:
            if controlDict[i]['curve_info']['type'] == 'PUMP1':
                print(i, 'Pump type 1...')

            elif controlDict[i]['curve_info']['type'] == 'PUMP2':
                print(i, 'Pump type 2...')

            elif controlDict[i]['curve_info']['type'] == 'PUMP3':
                head = h2 - h1 # pump pushes water from low head (h1) to higher head (h2)
                
                # x: Head
                # y: Flow
                x = [float(j) for j in controlDict[i]['curve_info']['x_val']] 
                y = [float(j) for j in controlDict[i]['curve_info']['y_val']]
                
                if head > max(x):
                    head = max(x)
                elif head < min(x):
                    head = min(x)

                # calculate q_full at given head.
                q_full = np.interp(
                    head,
                    np.array(x,dtype='float64'),
                    np.array(y,dtype='float64'),
                    )

                if controlDict[i]['q_goal'] == 0.0:
                    controlDict[i]['action'] = 0.0
                else:
                    controlDict[i]['action'] = controlDict[i]['q_goal'] / q_full

            elif controlDict[i]['curve_info']['type'] == 'PUMP4':
                print(i, 'Pump type 4...')
        
        
        # if target setting greater than 1, only open to 1.
        controlDict[i]['action'] = min(max(controlDict[i]['action'],0.0),1.0)
        controlDict[i]['pyswmmVar'].target_setting = controlDict[i]['action']
        controlDict[i]['actionTS'].append(controlDict[i]['action'])
        # controlDict['ORIFICE33@17312-17313']
        
        
def run_control_sim(control,controlDict,upstreamDict,downstreamDict,dsKeys,swmmINP,performanceDict,timestep):
    with Simulation(swmmINP) as sim:
        units = sim.system_units
        price = []
        PDemand = []
        
        nodes = Nodes(sim)
        links = Links(sim)
        
        
        # Initialize pyswmm variables for control, upstream, and downstream points.
        for point in controlDict:
            controlDict[point]['pyswmmVar'] = links[point]
        for point in upstreamDict:
            try:
                upstreamDict[point]['pyswmmVar'] = links[point]
            except:
                upstreamDict[point]['pyswmmVar'] = nodes[point]
        for point in downstreamDict:
            if downstreamDict[point]['type'] == 'link':
                downstreamDict[point]['pyswmmVar'] = links[point]
            if downstreamDict[point]['type'] == 'storage':
                downstreamDict[point]['pyswmmVar'] = nodes[point]
        
        # Initialize pyswmm variables for variables associated with performance metrics.
        for point in performanceDict:
            # IF node
            if performanceDict[point]['type'] == 'outfall' or performanceDict[point]['type'] == 'junction' or performanceDict[point]['type'] == 'storage':
                performanceDict[point]['pyswmmVar'] = nodes[point]
            # IF link
            elif performanceDict[point]['type'] == 'link' or performanceDict[point]['type'] == 'orifice':
                performanceDict[point]['pyswmmVar'] = links[point]
            else:
                pass
        
        print('running simulation...')
        for step in sim:
            
            
            # ---- Could make this subsection a function ---- 
            
            # append measures to timeseries for each upstream and downstream location
            for point in upstreamDict:
                upstreamDict[point]['ts'].append(upstreamDict[point]['pyswmmVar'].depth / upstreamDict[point]['max_depth'])
            for point in downstreamDict:
                try:
                    # if link
                    downstreamDict[point]['ts_flow'].append(downstreamDict[point]['pyswmmVar'].flow / downstreamDict[point]['max_flow'])
                    downstreamDict[point]['ts_depth'].append(downstreamDict[point]['pyswmmVar'].depth / downstreamDict[point]['max_depth'])
                except: 
                    # if node
                    downstreamDict[point]['ts_depth'].append(downstreamDict[point]['pyswmmVar'].depth / downstreamDict[point]['max_depth'])
            
            # append measures to timeseries for performance metric locations
            for point in performanceDict:
                # print(point)
                try: # if link
                    performanceDict[point]['ts_flow'].append(performanceDict[point]['pyswmmVar'].flow)
                except:  # if node
                    performanceDict[point]['ts_flow'].append(performanceDict[point]['pyswmmVar'].total_inflow)
            
            # ---- End subsection ----
            
            
            # try new function for mbc.
            if control:
                new_fnctn(upstreamDict, downstreamDict, controlDict, dsKeys, price, PDemand, nodes, timestep, units)
                
                # controlDict['ORIFICE33@17312-17313']['pyswmmVar'].target_setting = 0.0
                    
    print('done ...')
    winsound.Beep(440, 1000) # (Hz, millisecond)
    return price,PDemand


# def run_no_control_sim(controlDict,upstreamDict,downstreamDict,dsKeys,swmmINP,performanceDict):
#     with Simulation(swmmINP) as sim:
#         freq = '12s'
#         timesteps = pd.date_range(sim.start_time,sim.end_time,freq=freq)
#         price = []
#         PDemand = []
        
#         nodes = Nodes(sim)
#         links = Links(sim)
        
#         for point in controlDict:
#             controlDict[point]['pyswmmVar'] = links[point]
#         for point in upstreamDict:
#             upstreamDict[point]['pyswmmVar'] = links[point]
#         for point in downstreamDict:
#             downstreamDict[point]['pyswmmVar'] = links[point]
#         for point in performanceDict:
#             if performanceDict[point]['type'] == 'outfall' or performanceDict[point] == 'junction' or performanceDict[point] == 'storage':
#                 performanceDict[point]['pyswmmVar'] = nodes[point]
#             elif performanceDict[point]['type'] == 'link':
#                 performanceDict[point]['pyswmmVar'] = links[point]

#         print('running simulation...')
#         for step in sim:
#             for point in controlDict:
#                 controlDict[point]['pyswmmVar'].target_setting = 1.0
#             for point in upstreamDict:
#                 upstreamDict[point]['ts'].append(upstreamDict[point]['pyswmmVar'].depth / upstreamDict[point]['max_depth'])
#             for point in downstreamDict:
#                 downstreamDict[point]['ts_flow'].append(downstreamDict[point]['pyswmmVar'].flow / downstreamDict[point]['max_flow'])
#                 downstreamDict[point]['ts_depth'].append(downstreamDict[point]['pyswmmVar'].depth / downstreamDict[point]['max_depth'])
#             for point in performanceDict:
#                 try:
#                     performanceDict[point]['ts_flow'].append(performanceDict[point]['pyswmmVar'].flow)
#                 except:
#                     pass

#                 try:
#                     performanceDict[point]['ts_flow'].append(performanceDict[point]['pyswmmVar'].total_inflow)
#                 except:
#                     pass
#     print('done ...')

# def fnctn(ustream, dstream, setpts, uparam, dparam, n_tanks, action,max_flow,controlDict):
#     p = (sum(uparam*ustream) + sum(dparam*(dstream-setpts)))/(1 + n_tanks)
#     PD = np.zeros(n_tanks)
#     for i in range(0,n_tanks):
#         PD[i] = max(-p + uparam[i]*ustream[i],0)
#     PS = sum(PD)
    
#     for i in range(0,n_tanks):
#         if PS == 0:
#             Qi = 0
#         else:
#             # Qi = PD[i]/PS*setpts[0] # setpts[0] assumed to be downstream flow setpoint
#             Qi = sum(PD[0:i+1])/PS*setpts[0]*max_flow # setpts[0] assumed to be downstream flow setpoint
#         if ustream[i] == 0:
#             action[i] = 0.5
#         else:
#             h2i = Qi/(1.0*1*np.sqrt(2*9.81*ustream[i]))
# #             h2i = Qi/(0.61*1*np.sqrt(2*9.81*ustream[i]))
#             action[i] = max(min(h2i/2,1.0),0.0)
# #         if ustream[i] > 0.95:
# #             action[i] = 1.0

#     return p, PD, PS, action





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
            try:
                metadata = metadata + [i,saveDict['downstreamDict'][i]['epsilon'],saveDict['downstreamDict'][i]['gamma'],saveDict['downstreamDict'][i]['max_flow'],saveDict['downstreamDict'][i]['max_depth'],saveDict['downstreamDict'][i]['set_point'],saveDict['downstreamDict'][i]['set_derivative']]
            except:
                metadata = metadata + [i,saveDict['downstreamDict'][i]['epsilon'],saveDict['downstreamDict'][i]['gamma'],'no max flow',saveDict['downstreamDict'][i]['max_depth'],saveDict['downstreamDict'][i]['set_point'],saveDict['downstreamDict'][i]['set_derivative']]
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

def viz_shit(f_nocontrol,f_control,figname, save, control):
    # Unpack the pickles
    with open(f_control,'rb') as fC, open(f_nocontrol,'rb') as fNC:
        data_control = pickle.load(fC)
        data_nocontrol = pickle.load(fNC)

    fig,axarr = plt.subplots(2,2, figsize=(20,15))

    # Plot Control Points onto the figures
    leg = []
    for p1,p2 in zip(data_control['upstreamDict'],data_nocontrol['upstreamDict']):
        axarr[0,0].plot(data_control['upstreamDict'][p1]['ts'])
        axarr[0,0].plot(data_nocontrol['upstreamDict'][p2]['ts'],linestyle='--',linewidth=3)
        leg = leg + [data_control['upstreamDict'][p1]['location'],data_nocontrol['upstreamDict'][p1]['location']]
    for p1,p2 in zip(data_control['downstreamDict'],data_nocontrol['downstreamDict']):
        axarr[0,0].plot(data_control['downstreamDict'][p1]['ts_depth'])
        axarr[0,0].plot(data_nocontrol['downstreamDict'][p2]['ts_depth'],color='k',linestyle='--')
        leg = leg + [data_control['downstreamDict'][p1]['location'],data_nocontrol['downstreamDict'][p1]['location']]
    axarr[0,0].set_title('Normalized Depth at Points of Interest')
    
    # make legend
    axarr[0,0].legend(leg,title = 'Solid = Control, Dashed = No Control', ncol=2)
    axarr[0,0].set_ylabel('Normalized Depth')
    
    leg = []
    t = np.array([i for i in range(0,len(data_control['downstreamDict'][p1]['ts_depth']))])
    if control:
        for i in data_control['controlDict']:
            axarr[1,0].scatter(t,data_control['controlDict'][i]['actionTS'])
            leg = leg + [data_control['controlDict'][i]['location']]
        axarr[1,0].set_ylabel('Action Setting')
        axarr[1,0].set_title('Actions at Control Points')
        axarr[1,0].legend(leg)
    
    
    axarr[0,1].plot(data_control['performanceDict']['1002']['ts_flow'])
    axarr[0,1].plot(data_nocontrol['performanceDict']['1002']['ts_flow'],color='k',linestyle='--',linewidth=3)
    axarr[0,1].set_title('WWTP Inflow')
    axarr[0,1].set_ylabel('Flow [cfs]')
    
    axarr[0,1].plot(data_control['performanceDict']['2908']['ts_flow'])
    axarr[0,1].plot(data_nocontrol['performanceDict']['2908']['ts_flow'],linestyle='--',linewidth=3)
    leg = ['WWTP Inflow','WWTP Inflow','CC RTB Overflow','CC RTB Overflow']
    axarr[0,1].legend(leg,title = 'Solid = Control, Dashed = No Control', ncol=2)
    axarr[0,1].set_title('Inlows and Overflows')
    
    
    
    leg = []
    for i in data_control['upstreamDict']:
        axarr[1,1].plot(data_control['upstreamDict'][i]['PD'])
        leg = leg + [data_control['upstreamDict'][i]['location']]
    axarr[1,1].plot(data_control['price'],color = 'k', linewidth = 3)
    leg = leg + ['Price']
    axarr[1,1].legend(leg,title = 'Downstream Price, Upstream Wealth')
    axarr[1,1].set_ylabel('Currency')
    axarr[1,1].set_title('The Market')
    
    
    if save:
        plt.savefig(figname)
    
    plt.show()
    plt.close()
    
    
    
    