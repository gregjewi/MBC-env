from collections import OrderedDict
import swmmAPImini as sm
import mbc

swmmINP = 'GDRSS_25YR24HR.inp'
simNotes = ['qa qc orifice equations']

pickleOut = 'orifice_qaqc.p'
metaCSV = 'CC_FV_FREUD_meta.csv'

cF = 'input_files/ControlPoints_Init_ISDs.csv'
uF = 'input_files/Upstream_Init.csv'
dF = 'input_files/Downstream_Init.csv'
pF = 'input_files/PerformancePoints_Init.csv' # Elements of interest in SWMM to analyze performance

controlDict = sm.return_inputs(cF,'control')
upstreamDict = sm.return_inputs(uF,'upstream')
downstreamDict = sm.return_inputs(dF,'downstream')
performanceDict = sm.return_inputs(pF,'performance')
dsKeys = list(downstreamDict.keys())

# With or without control?
control = True

conduits,nodes,storages,subcatchments,outfalls,orifices,pumps = sm.make_element_dictionaries()

sm.get_depth(upstreamDict,conduits,storages)
sm.get_depth(downstreamDict,conduits,storages)
sm.get_q_full(downstreamDict,conduits)
sm.orifice_xsect_grab(controlDict,orifices)
sm.pump_curve_grab(controlDict,pumps)

sm.performance_elements(performanceDict,conduits,nodes,storages,subcatchments,outfalls)
dsKeys = list(downstreamDict.keys())

try:
    del conduits
    del nodes
    del storages
    del subcatchments
    del outfalls
except:
    pass

downstreamDict[dsKeys[0]]['max_flow'] = 1.0

if control:
    timesteps, price, PDemand = mbc.run_control_sim(controlDict,upstreamDict,downstreamDict,dsKeys,swmmINP,performanceDict)
else:
    mbc.run_no_control_sim(controlDict,upstreamDict,downstreamDict,dsKeys,swmmINP,performanceDict)


saveDict = {
    'controlDict':controlDict,
    'upstreamDict':upstreamDict,
    'downstreamDict':downstreamDict,
    'metaCSV':metaCSV,
    'pickleOut':pickleOut,
    'notes':simNotes,
    'file':swmmINP,
    'control':str(control),
    'performanceDict':performanceDict
}

try:
    saveDict['price'] = price
    saveDict['PDemand'] = PDemand
except:
    pass

mbc.save_that_shit(saveDict)