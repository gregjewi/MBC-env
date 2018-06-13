from collections import OrderedDict
import swmmAPImini as sm
import mbc

swmmINP = 'NoWestside_withISDs_RainEvents.inp'


# With or without control?
control = True

cF = 'input_files/ControlPoints_Init_CC_FV_FREUD.csv'
uF = 'input_files/Upstream_Init_CC_FV_FREUD.csv'
dF = 'input_files/Downstream_Init_CC_FV_FREUD.csv'
pF = 'input_files/PerformancePoints_Init.csv' # Elements of interest in SWMM to analyze performance

conduits,nodes,storages,subcatchments,outfalls,orifices,pumps = sm.make_element_dictionaries()

metaCSV = 'CC_FV_FREUD_meta.csv'

for a in range(4,10):
	for b in range(0,10):
		pickleOut = 'one_inch_' + str(a) + '_' + str(b) + '.p'
		simNotes = ['for loops a:' + str(a) + ' b: ' + str(b)]

		controlDict = sm.return_inputs(cF,'control')
		upstreamDict = sm.return_inputs(uF,'upstream')
		downstreamDict = sm.return_inputs(dF,'downstream')
		performanceDict = sm.return_inputs(pF,'performance')
		dsKeys = list(downstreamDict.keys())

		sm.get_depth(upstreamDict,conduits,storages)
		sm.get_depth(downstreamDict,conduits,storages)
		sm.get_q_full_and_other(downstreamDict,conduits,storages)
		sm.orifice_xsect_grab(controlDict,orifices)
		sm.pump_curve_grab(controlDict,pumps)

		timestep = sm.get_timestep(swmmINP)

		sm.performance_elements(performanceDict,conduits,nodes,storages,subcatchments,outfalls,orifices)
		dsKeys = list(downstreamDict.keys())

		upstreamDict['17311']['uparam'] = a
		upstreamDict['5220']['uparam'] = a
		upstreamDict['5010']['uparam'] = a
		upstreamDict['5010']['max_depth'] = 24.0
		downstreamDict['2909']['max_depth'] = 22.0
		downstreamDict['2909']['set_point'] = 0.8
		downstreamDict['2909']['epsilon'] = b
		downstreamDict['2909']['gamma'] = 0.0

		price, PDemand = mbc.run_control_sim(control,controlDict,upstreamDict,downstreamDict,dsKeys,swmmINP,performanceDict,timestep)

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

		del upstreamDict
		del controlDict
		del downstreamDict
		del price
		del PDemand
		del performanceDict