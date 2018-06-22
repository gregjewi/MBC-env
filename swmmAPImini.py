import math
import pandas as pd
from collections import OrderedDict
from pyswmm import Simulation

def make_element_dictionaries():
    conduitF = 'swmmModelElements/CONDUITS.txt'
    xsF = 'swmmModelElements/XSECTIONS.txt'
    jF = 'swmmModelElements/JUNCTIONS.txt'
    storF = 'swmmModelElements/STORAGE.txt'
    subF = 'swmmModelElements/SUBCATCHMENTS.txt'
    outfallF = 'swmmModelElements/OUTFALLS.txt'
    orificeF = 'swmmModelElements/ORIFICES.txt'
    pumpF = 'swmmModelElements/PUMPS.txt'
    curvesF = 'swmmModelElements/CURVES.txt'
    
    conduits = make_conduit_dictionary(conduitF,xsF)
    nodes = make_node_dictionary(jF)
    storages = make_storage_dictionary(storF)
    conduits = calc_slope(conduits,nodes,storages,0)
    conduits = calc_qfull(conduits)
    subcatchments = make_subcatchment_dictionary(subF)
    outfalls = make_outfall_dictionary(outfallF)
    orifices = make_orifice_dictionary(orificeF,xsF)
    curves = make_curves_dictionary(curvesF)
    pumps = make_pump_dictionary(pumpF,curves)
    datum_conversion(nodes,479.755) #pass dict and offset
    datum_conversion(storages,479.755) #pass dict and offset
    datum_conversion(outfalls,479.755) #pass dict and offset

    system = {}

    return conduits,nodes,storages,subcatchments,outfalls,orifices,pumps

def make_conduit_dictionary(conduitF,xsF):
    conduitDict = {}
    myList = []
    
    with open(conduitF) as f:
        for l in f:
            if l[0].isalnum():
                myList.append(l.split())
    
    for l in myList:
        conduitDict[l[0]] = {
            'from_node': l[1],
            'to_node': l[2],
            'length': float(l[3]),
            'roughness': float(l[4]),
            'in_offset': float(l[5]),
            'out_offset': float(l[6]),
            'init_flow': float(l[7]),
            'max_flow': float(l[8]),
        }
    
    myList = []
    with open(xsF) as f:
        for l in f:
            if l[0].isalnum():
                myList.append(l.split())
    
    for l in myList:
        try:
            conduitDict[l[0]]['shape'] = l[1]
            conduitDict[l[0]]['geom1'] = float(l[2])
            conduitDict[l[0]]['geom2'] = float(l[3])
            conduitDict[l[0]]['geom3'] = float(l[4])
            conduitDict[l[0]]['geom4'] = float(l[5])
            conduitDict[l[0]]['barrels'] = int(l[6])
            conduitDict[l[0]]['culvert'] = float(l[7])
        except:
            pass

    return conduitDict

def make_node_dictionary(jF):
    nodeDict = {}
    myList = []

    with open(jF) as f:
        for l in f:
            if l[0].isalnum():
                myList.append(l.split())
                
    for l in myList:
        nodeDict[l[0]] = {
            'elevation': float(l[1]),
            'max_depth': float(l[2]),
            'init_depth': float(l[3]),
            'sur_depth': float(l[4]),
            'a_ponded': float(l[5])
        }
    
    return nodeDict

def make_storage_dictionary(storF):
    storageDict = {}
    myList = []
    
    with open(storF) as f:
        for l in f:
            if l[0].isalnum():
                myList.append(l.split())
    
    for l in myList:
        storageDict[l[0]] = {
            'elevation':float(l[1]),
            'max_depth':float(l[2]),
            'init_depth:':float(l[3]),
            'shape':l[4],
            #Curve Name/Params
            #N/A
            #Fevap
            #PSI
            #Ksat
            #IMD
        }
        
        if l[4] == 'FUNCTIONAL':
            storageDict[l[0]]['A'] = float(l[5])
            storageDict[l[0]]['B'] = float(l[6])
            storageDict[l[0]]['C'] = float(l[7])
            
    storageDict = calc_storage_vol(storageDict)
    
    return storageDict

def make_subcatchment_dictionary(subF):
    myList = []
    Subcatchment = {}
    with open(subF,'r') as f:
        for l in f:
            if l[0].isalnum():
                myList.append(l.split())

        for l in myList:
            Subcatchment[l[0]] = {
                'rain_gage': l[1],
                'outlet': l[2],
                'area': float(l[3]),
                'per_imperv': float(l[4]),
                'width': float(l[5]),
                'slope': float(l[6]),
                'curblen': float(l[7])
            }
            
    return Subcatchment

def make_outfall_dictionary(outfallF):
    myList = []
    outfall = {}
    
    with open(outfallF,'r') as f:
        for l in f:
            
            if l[0].isalnum():
                myList.append(l.split())
                
        for l in myList:
            outfall[l[0]] = {
                'elevation': float(l[1]),
                'type': l[2],
                'stage_data': float(l[3]),
                'gated': l[4],
                
        # try later "route to"
            }
            
    return outfall

def make_orifice_dictionary(orificeF,xsF):
    myList = []
    orifice = {}
    
    with open(orificeF,'r') as f:
        for l in f:
            if l[0].isalnum():
                myList.append(l.split())
        
        for l in myList:
            orifice[l[0]] = {
                'from_node': str(l[1]),
                'to_node': str(l[2]),
                'type':l[3],
                'offset':float(l[4]),
                'Cd':float(l[5]),
                'gated':l[6],
                'close_time':float(l[7])
            }
            
    myList = []
    with open(xsF) as f:
        for l in f:
            if l[0].isalnum():
                myList.append(l.split())
    
    for l in myList:
        try:
            orifice[l[0]]['shape'] = l[1]
            orifice[l[0]]['geom1'] = float(l[2])
            orifice[l[0]]['geom2'] = float(l[3])
            orifice[l[0]]['geom3'] = float(l[4])
            orifice[l[0]]['geom4'] = float(l[5])
            orifice[l[0]]['barrels'] = int(l[6])
            orifice[l[0]]['culvert'] = float(l[7])
        except:
            pass
    
    return orifice

def make_pump_dictionary(pumpF,curves):
    myList = []
    pump = {}
    
    with open(pumpF,'r') as f:
        for l in f:
            if l[0].isalnum():
                myList.append(l.split())
            
        for l in myList:
            pump[l[0]] = {
                'from_node':l[1],
                'to_node':l[2],
                'pump_curve':l[3],
                'status':l[4],
                'startup':l[5],
                'shutoff':l[6]
            }
    
    for p in pump:
        pump[p]['curve_info'] = curves[pump[p]['pump_curve']]
    
    return pump

def make_curves_dictionary(curvesF):
    myList = []
    curves = {}
    
    with open(curvesF,'r') as f:
        for l in f:
            if l[0].isalnum():
                myList.append(l.split())
                
    for l in myList:
        if len(l) == 4:
            curves[l[0]] = {
                'type': l[1],
                'x_val': [l[2]],
                'y_val': [l[3]]
            }
        if len(l) == 3:
            curves[l[0]]['x_val'].append(l[1])
            curves[l[0]]['y_val'].append(l[2])
            
    return curves

def calc_storage_vol(storageDict):
    for element in storageDict:
        if storageDict[element]['shape'] == 'FUNCTIONAL':
            storageDict[element]['total_storage'] = ( storageDict[element]['A'] * storageDict[element]['max_depth']**(storageDict[element]['B'] + 1) / ( storageDict[element]['B'] + 1 ) ) + ( storageDict[element]['C'] * storageDict[element]['max_depth'] )
        else:
            # do something here to add in storage curves
            pass
    
    return storageDict

def calc_slope(conduitDict,nodeDict,storageDict,min_slope):
    for item in conduitDict:
        if conduitDict[item]['from_node'] in nodeDict.keys():
            e1 = nodeDict[conduitDict[item]['from_node']]['elevation']+conduitDict[item]['in_offset']
        elif conduitDict[item]['from_node'] in storageDict.keys():
            e1 = storageDict[conduitDict[item]['from_node']]['elevation']+conduitDict[item]['in_offset']
        else:
            e1 = 1
        
        if conduitDict[item]['to_node'] in nodeDict.keys():
            e2 = nodeDict[conduitDict[item]['to_node']]['elevation']+conduitDict[item]['out_offset']
        elif conduitDict[item]['to_node'] in storageDict.keys():
            e2 = storageDict[conduitDict[item]['to_node']]['elevation']+conduitDict[item]['out_offset']
        else:
            e2 = 1
            
        if e1==1 or e2==1:
            conduitDict[item]['slope_flag'] = True
        else:
            conduitDict[item]['slope_flag'] = False
        
        slope = (e1 - e2)/conduitDict[item]['length']
        
        if slope < min_slope:
            slope = min_slope
            conduitDict[item]['slope_flag'] = True
        
        conduitDict[item]['slope'] = slope
    
    return conduitDict

def calc_qfull(conduits):
    for item in conduits:
        if conduits[item]['shape'] == 'CIRCULAR':
            # compute Qfull as pipe full manning equation
            conduits[item]['q_full'] = (conduits[item]['geom1']**(8/3)*conduits[item]['slope']**(1/2))/(4**(5/3)*conduits[item]['roughness'])*math.pi
        elif conduits[item]['shape'] == 'RECT_CLOSED':
            # Compute q_full as manning equation of pipe with manning eq with depth as 0.95
            conduits[item]['q_full'] = (1.49/conduits[item]['roughness']) * (0.95 * conduits[item]['geom1'] * conduits[item]['geom2']) * (conduits[item]['geom2'] * 0.95 * conduits[item]['geom1'] / (conduits[item]['geom2'] + 2 * 0.95 * conduits[item]['geom1']))**(2/3)
        else:
            conduits[item]['q_full'] = 1;
            
    return conduits

def change_elev(swmm_nodes,nodeD,storD,outD):
    for n in swmm_nodes:
        try:
            #print(swmm_nodes[n.nodeid].invert_elevation)
            swmm_nodes[n.nodeid].invert_elevation = nodeD[n.nodeid]['elev_detroit']
            #print(swmm_nodes[n.nodeid].invert_elevation)
        except:
            try:
                swmm_nodes[n.nodeid].invert_elevation = storD[n.nodeid]['elev_detroit']
                #print('beech')
            except:
                try:
                    swmm_nodes[n.nodeid].invert_elevation = outD[nodeid]['elev_detroit']
                    #print('grab em')
                except:
                    print("Error in Datum Change for "+n.nodeid)
                    pass
        
def datum_conversion(withElev,offset):
    for point in withElev:
        withElev[point]['elev_detroit'] = withElev[point]['elevation'] - offset

def return_inputs(filename,handle):
    #
    df = pd.read_csv(filename,index_col='name')
    df.index = df.index.map(str)
    df_dict = df.to_dict(into=OrderedDict,orient='index')

    for i in df_dict:
        if handle == 'control':
            df_dict[i]['actionTS'] = []
            df_dict[i]['q_goal'] = 1.0
        elif handle == 'upstream':
            df_dict[i]['ts'] = []
            df_dict[i]['PD'] = []
        elif handle == 'downstream':
            df_dict[i]['ts_flow'] = []
            df_dict[i]['ts_depth'] = []
        elif handle == 'performance':
            df_dict[i]['ts_flow'] = []
            df_dict[i]['ts_depth'] = []
        else:
            pass

    return df_dict

def get_depth(elements,conduits,storages):
    for element in elements:
        if elements[element]['type'] == 'link':
            elements[element]['max_depth'] = conduits[element]['geom1']
        elif elements[element]['type'] == 'storage':
            elements[element]['max_depth'] = storages[element]['max_depth']
        else:
            pass

def get_q_full_and_other(elements,conduits,storages):
    for element in elements:
        if elements[element]['type'] == 'link':
            elements[element]['max_flow'] = conduits[element]['q_full']
        elif elements[element]['type'] == 'storage':
            elements[element]['total_storage'] = storages[element]['total_storage']
        else:
            pass

def performance_elements(elements,conduits,nodes,storages,subcatchments,outfalls,orifices):
    for element in elements:
        if elements[element]['type'] == 'outfall':
            elements[element]['elevation'] = outfalls[element]['elevation']
            elements[element]['elev_detroit'] = outfalls[element]['elev_detroit']
        elif elements[element]['type'] == 'link':
            elements[element]['max_depth'] = conduits[element]['geom1']
        elif elements[element]['type'] == 'storage':
            elements[element]['max_depth'] = storages[element]['max_depth']
        elif elements[element]['type'] == 'orifice':
            elements[element]['max_depth'] = orifices[element]['geom1']
        else:
            pass

def orifice_xsect_grab(controlDict,orifices):
    # Add items from orifices dictionary to the control dict
    for i in controlDict:
        try:
            controlDict[i].update(orifices[i])
        except:
            pass
        
def pump_curve_grab(controlDict, pumps):
    # Add information to controlDict with pump curve info.
    for i in controlDict:
        try:
            controlDict[i].update(pumps[i])
        except:
            pass
        
def get_timestep(inF):
    with open(inF) as inputF:
        routing_step = False
        for l in inputF:
            if l[:12] == 'ROUTING_STEP':
                routing_step = l
            else:
                pass

            if routing_step:
                step_str = routing_step.split()[1].split(':')
                
                timestep_sec = float(step_str[0])*60*60 + float(step_str[1])*60 + float(step_str[2])
                
                return timestep_sec
    