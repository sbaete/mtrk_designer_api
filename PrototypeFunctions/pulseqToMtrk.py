################################################################################
### mtrk project - Pypulseq-based conversion tool from SDL to Pulseq.        ###
### Version 0.1.1                                                            ###
### Anais Artiges and the mtrk project team at NYU - 04/29/2024              ###
################################################################################

import numpy as np
from numpy import linspace
from scipy import interpolate
from scipy.integrate import simpson
import sys
import math 
import pypulseq
from pypulseq.event_lib import EventLibrary
from pypulseq.supported_labels_rf_use import get_supported_labels
from types import SimpleNamespace
import matplotlib.pyplot as plt

from SDL_read_write.pydanticSDLHandler import *

def pulseqToMtrk(file_path):
    sequence_data = PulseSequence(file = File(), 
                                  settings = Settings(), 
                                  infos = Info(),
                                  instructions = {}, 
                                  objects = {}, 
                                  arrays={}, 
                                  equations={})
    seq = pypulseq.Sequence()
    readPulseq(seq, file_path)
    # print("definitions", seq.definitions)
    # print("blockEvents", seq.blockEvents)
    # print("gradLibrary", seq.gradLibrary)
    # print("rfLibrary", seq.rfLibrary)
    # print("adcLibrary", seq.adcLibrary)
    # print("trigLibrary", seq.trigLibrary)
    # print("labelsetLibrary", seq.labelsetLibrary)
    # print("labelincLibrary", seq.labelincLibrary)
    # print("extensionStringIDs", seq.extensionStringIDs)
    # print("extensionNumericIDs", seq.extensionNumericIDs)
    # print("blockDurations", seq.blockDurations)
    # print("gradRasterTime", seq.gradRasterTime)
    # TO DO correct the mtrkToPulseq for correct raster time
    seq.rfRasterTime = 20e-6
    print("rfRasterTime", seq.rfRasterTime)
    # print("adcRasterTime", seq.adcRasterTime)
    # print("blockDurationRaster", seq.blockDurationRaster)
    # print("signatureType", seq.signatureType)
    # print("signatureValue", seq.signatureValue)
    # print("signatureFile", seq.signatureFile)
    # print("shapeLibrary", seq.shapeLibrary)

    # for counter in range(0, len(seq.blockEvents)):
    # rf_events = []
    # gx_events = []
    # gy_events = []
    # gz_events = []
    # adc_events = []
    # # print("seq.blockEvents", seq.blockEvents)
    # for event in seq.blockEvents:
    #     # print("event", event)
    #     rf_events.append(event[1])
    #     gx_events.append(event[2])
    #     gy_events.append(event[3])
    #     gz_events.append(event[4])
    #     adc_events.append(event[5])
    # print("rf_events", rf_events)

    periodicEvents, variableEvents, variableIndexes = evaluatePeriodicity(seq.blockEvents)
    extractRFevents(seq, periodicEvents, variableEvents, variableIndexes)
    return seq

def readPulseq(obj, filename, *args):
    # READ Load sequence from file.
    #   READ(seqObj, filename, ...) Read the given filename and load sequence
    #   data into sequence object.
    #
    #   optional parameter 'detectRFuse' can be given to let the function
    #   infer the currently missing flags concerning the intended use of the RF
    #   pulses (excitation, refocusing, etc). These are important for the
    #   k-space trajectory calculation
    #
    #   Examples:
    #   Load the sequence defined in gre.seq in my_sequences directory
    #
    #       read(seqObj,'my_sequences/gre.seq')
    #
    # See also  write

    detectRFuse = False
    if args and 'detectRFuse' in args:
        detectRFuse = True

    fid = open(filename)

    # Clear sequence data
    obj.blockEvents = []
    obj.definitions = {}
    obj.gradLibrary = EventLibrary()
    obj.shapeLibrary = EventLibrary()
    obj.rfLibrary = EventLibrary()
    obj.adcLibrary = EventLibrary()
    obj.trigLibrary = EventLibrary()
    obj.labelsetLibrary = EventLibrary()
    obj.labelincLibrary = EventLibrary()
    obj.extensionStringIDs = []
    obj.extensionNumericIDs = []

    version_combined = 0

    # Load data from file
    while True:
        section = ""
        lineData = fid.readline().strip()
        if lineData != "" and lineData[0] == '[':
            section = lineData

        if section == '[DEFINITIONS]':
            obj.definitions = readDefinitions(fid)
            v = obj.get_definition('GradientRasterTime')
            if v:
                obj.gradRasterTime = v
            v = obj.get_definition('RadiofrequencyRasterTime')
            if v:
                obj.rfRasterTime = v
            v = obj.get_definition('AdcRasterTime')
            if v:
                obj.adcRasterTime = v
            v = obj.get_definition('BlockDurationRaster')
            if v:
                obj.blockDurationRaster = v
        elif section == '[SIGNATURE]':
            tmpSignDefs = readDefinitions(fid)
            if 'Type' in tmpSignDefs:
                obj.signatureType = tmpSignDefs['Type']
            if 'Hash' in tmpSignDefs:
                obj.signatureValue = tmpSignDefs['Hash']
                obj.signatureFile = 'Text'  # we are reading a text file, so much is known for sure
                break
        elif section == '[VERSION]':
            version_major, version_minor, version_revision = readVersion(fid)
            assert version_major == obj.version_major, 'Unsupported version_major {}'.format(version_major)
            version_combined = 1000000 * version_major + 1000 * version_minor + version_revision
            if version_combined < 1002000:
                raise ValueError('Unsupported version {}.{}.{}, only file format revision 1.2.0 and above are supported'.format(version_major, version_minor, version_revision))
            if version_combined < 1003001:
                print('Loading older Pulseq format file (version {}.{}.{}), some code may function not as expected'.format(version_major, version_minor, version_revision))
        elif section == '[BLOCKS]':
            if 'version_major' not in locals():
                raise ValueError('Pulseq file MUST include [VERSION] section prior to [BLOCKS] section')
            obj.blockEvents, obj.blockDurations, delayInd_tmp = readBlocks(fid, obj.blockDurationRaster, version_combined)
        elif section == '[RF]':
            if version_combined >= 1004000:
                obj.rfLibrary = readEvents(fid, [1, 1, 1, 1, 1e-6, 1, 1])  # this is 1.4.x format
            else:
                obj.rfLibrary = readEvents(fid, [1, 1, 1, 1e-6, 1, 1])  # this is 1.3.x and below
        elif section == '[GRADIENTS]':
            if version_combined >= 1004000:
                obj.gradLibrary = readEvents(fid, [1, 1, 1, 1e-6], 'g', obj.gradLibrary)  # this is 1.4.x format
            else:
                obj.gradLibrary = readEvents(fid, [1, 1, 1e-6], 'g', obj.gradLibrary)  # this is 1.3.x and below
        elif section == '[TRAP]':
            obj.gradLibrary = readEvents(fid, [1, 1e-6, 1e-6, 1e-6, 1e-6], 't', obj.gradLibrary)
        elif section == '[ADC]':
            obj.adcLibrary = readEvents(fid, [1, 1e-9, 1e-6, 1, 1])
        elif section == '[DELAYS]':
            if version_combined >= 1004000:
                raise ValueError('Pulseq file revision 1.4.0 and above MUST NOT contain [DELAYS] section')
            tmp_delayLibrary = readEvents(fid, 1e-6)
        elif section == '[SHAPES]':
            obj.shapeLibrary = readShapes(fid, False)
        elif section == '[EXTENSIONS]':
            obj.extensionLibrary = readEvents(fid)
        elif section == "":
            pass
        else:
            if section.startswith('extension TRIGGERS'):
                id = int(section[19:])
                obj.setExtensionStringAndID('TRIGGERS', id)
                obj.trigLibrary = readEvents(fid, [1, 1, 1e-6, 1e-6])
            elif section.startswith('extension LABELSET'):
                id = int(section[19:])
                obj.setExtensionStringAndID('LABELSET', id)
                obj.labelsetLibrary = readAndParseEvents(fid, int, get_supported_labels())
            elif section.startswith('extension LABELINC'):
                id = int(section[19:])
                obj.setExtensionStringAndID('LABELINC', id)
                obj.labelincLibrary = readAndParseEvents(fid, int, get_supported_labels())
            else:
                raise ValueError('Unknown section code: {}'.format(section))

    fid.close()

    # if version_combined < 1002000:
    #     raise ValueError('Unsupported version {:07d}, only file format revision 1.2.0 (1002000) and above are supported'.format(version_combined))

    # if version_combined < 1004000:
    #     # fix blocks, gradients and RF objects imported from older versions
    #     for i in range(len(obj.rfLibrary.data)):
    #         obj.rfLibrary.data[i].array = obj.rfLibrary.data[i].array[:3] + [0] + obj.rfLibrary.data[i].array[3:]
    #         obj.rfLibrary.lengths[i] += 1

    #     for i in range(len(obj.gradLibrary.data)):
    #         if obj.gradLibrary.type[i] == 't':
    #             if obj.gradLibrary.data[i].array[2] == 0:
    #                 if abs(obj.gradLibrary.data[i].array[1]) == 0 and obj.gradLibrary.data[i].array[3] > 0:
    #                     obj.gradLibrary.data[i].array[3] -= obj.gradRasterTime
    #                     obj.gradLibrary.data[i].array[2] = obj.gradRasterTime
    #             if obj.gradLibrary.data[i].array[4] == 0:
    #                 if abs(obj.gradLibrary.data[i].array[1]) == 0 and obj.gradLibrary.data[i].array[3] > 0:
    #                     obj.gradLibrary.data[i].array[3] -= obj.gradRasterTime
    #                     obj.gradLibrary.data[i].array[4] = obj.gradRasterTime
    #         if obj.gradLibrary.type[i] == 'g':
    #             obj.gradLibrary.data[i].array = obj.gradLibrary.data[i].array[:2] + [0] + obj.gradLibrary.data[i].array[3:]
    #             obj.gradLibrary.lengths[i] += 1

    #     obj.blockDurations = [0] * len(obj.blockEvents)
    #     for iB in range(len(obj.blockEvents)):
    #         b = obj.get_block(iB)
    #         if delayInd_tmp[iB] > 0:
    #             b.delay.type = 'delay'
    #             b.delay.delay = tmp_delayLibrary.data[delayInd_tmp[iB]].array
    #         obj.blockDurations[iB] = pypulseq.calc_duration(b)

    # gradChannels = ['gx', 'gy', 'gz']
    # gradPrevLast = [0] * len(gradChannels)
    # for iB in range(len(obj.blockEvents)):
    #     print("obj.blockEvents[iB]",obj.blockEvents[iB])
    #     print("iB",iB)
    #     obj.add_block(obj.blockEvents[iB])
    #     b = obj.get_block(iB)
    #     block_duration = obj.blockDurations[iB]
    #     eventIDs = obj.blockEvents[iB]
    #     for j in range(len(gradChannels)):
    #         grad = b[gradChannels[j]]
    #         if not grad:
    #             gradPrevLast[j] = 0
    #             continue
    #         if grad.type == 'grad':
    #             if grad.delay > 0:
    #                 gradPrevLast[j] = 0
    #             if hasattr(grad, 'first'):
    #                 continue
    #             grad.first = gradPrevLast[j]
    #             if grad.time_id != 0:
    #                 grad.last = grad.waveform[-1]
    #                 grad_duration = grad.delay + grad.tt[-1]
    #             else:
    #                 odd_step1 = [grad.first] + [2 * x for x in grad.waveform]
    #                 odd_step2 = [x * ((i % 2) * 2 - 1) for i, x in enumerate(odd_step1)]
    #                 waveform_odd_rest = [x * ((i % 2) * 2 - 1) for i, x in enumerate(cumsum(odd_step2))]
    #                 grad.last = waveform_odd_rest[-1]
    #                 grad_duration = grad.delay + len(grad.waveform) * obj.gradRasterTime
    #             gradPrevLast[j] = grad.last
    #             if grad_duration + eps < block_duration:
    #                 gradPrevLast[j] = 0
    #             id = eventIDs[j + 2]
    #             amplitude = obj.gradLibrary.data[id].array[0]
    #             if version_combined >= 1004000:
    #                 old_data = [amplitude, grad.shape_id, grad.time_id, grad.delay]
    #             else:
    #                 old_data = [amplitude, grad.shape_id, grad.delay]
    #             new_data = [amplitude, grad.shape_id, grad.time_id, grad.delay, grad.first, grad.last]
    #             update_data(obj.gradLibrary, id, old_data, new_data, 'g')
    #         else:
    #             gradPrevLast[j] = 0

    if detectRFuse:
        for k in obj.rfLibrary.keys:
            libData = obj.rfLibrary.data[k].array
            rf = obj.rfFromLibData(libData)
            flipAngleDeg = abs(sum(rf.signal[:-1] * (rf.t[1:] - rf.t[:-1]))) * 360
            offresonance_ppm = 1e6 * rf.freqOffset / obj.sys.B0 / obj.sys.gamma
            if flipAngleDeg < 90.01:
                obj.rfLibrary.type[k] = 'e'
            else:
                if rf.shape_dur > 6e-3 and -3.5 <= offresonance_ppm <= -3.4:
                    obj.rfLibrary.type[k] = 's'
                else:
                    obj.rfLibrary.type[k] = 'r'

    return

def readDefinitions(fid):
    def_dict = {}
    line = fid.readline().strip()
    while line and not (line.isspace()):
        if line.startswith('#'):
            line = fid.readline().strip()
        else:
            tok = line.split()
            tok[1] = tok[1].replace("e", "E")
            key = tok[0]
            if is_float(tok[1]) == False:  
                 value = tok[1]
            else:
                value = float(tok[1]) if len(tok) > 1 else line[len(key) + 1:].strip()
            def_dict[key] = value
            line = fid.readline().strip()
    return def_dict

def readVersion(fid):
    major = None
    minor = None
    revision = None
    line = fid.readline().strip()
    while line and not (line.isspace() or line.startswith('#')):
        tok = line.split()
        if tok[0] == 'major':
            major = int(tok[1])
        elif tok[0] == 'minor':
            minor = int(tok[1])
        elif tok[0] == 'revision':
            revision = int(tok[1])
        line = fid.readline().strip()
    return major, minor, revision

def readBlocks(fid, blockDurationRaster, version_combined):
    eventTable = []
    blockDurations = []
    delayIDs_tmp = []
    line = fid.readline().strip()
    ### Warning : the .append() might be an issue as the original code was adding at index blockEvent[0] instead
    while line and not (line.isspace() or line.startswith('#')):
        blockEvents = [int(i) for i in line.split()]
        eventTable.append([0] + blockEvents[2:] + [0]) if version_combined <= 1002001 else eventTable.append([0] + blockEvents[2:])
        if version_combined >= 1004000:
            blockDurations.append(blockEvents[1] * blockDurationRaster)
        else:
            delayIDs_tmp.append(blockEvents[1])
        line = fid.readline().strip()
    return eventTable, blockDurations, delayIDs_tmp

def readEvents(fid, scale, type=None, eventLibrary=None):
    if not scale:
        scale = [1]
    if not eventLibrary:
        eventLibrary = EventLibrary()
    line = fid.readline().strip()
    while line and not (line.isspace() or line.startswith('#')):
        data = list(map(float, line.split()))
        id = int(data[0])
        data = [x * y for x, y in zip(scale, data[1:])]
        if type:
            eventLibrary.insert(id, data, type)
        else:
            eventLibrary.insert(id, data)
        line = fid.readline().strip()
    return eventLibrary

def readAndParseEvents(fid, *parsers):
    eventLibrary = EventLibrary()
    line = fid.readline().strip()
    while line and not (line.isspace() or line.startswith('#')):
        datas = line.split()
        data = []
        id = int(datas[0])
        for i in range(1, len(datas)):
            if i > len(parsers):
                data.append(int(datas[i]))
            else:
                data.append(parsers[i - 1](datas[i]))
        eventLibrary.insert(id, data)
        line = fid.readline().strip()
    return eventLibrary

def readShapes(fid, forceConvertUncompressed):
    shapeLibrary = EventLibrary()
    line = fid.readline().strip()
    line = fid.readline().strip()
    while line and not line.isspace():
        tok = line.split()
        id = int(tok[1])
        line = fid.readline().strip()
        tok = line.split()
        num_samples = int(tok[1])
        data = []
        line = fid.readline().strip()
        while line and not (line.isspace() or line.startswith('#')):
            data.extend(list(map(float, line.split())))
            line = fid.readline().strip()

        line = fid.readline().strip()

        if len(data) != num_samples:
            compressed_shape = SimpleNamespace(num_samples = num_samples, data = np.array(data))
            shape = list(pypulseq.decompress_shape.decompress_shape(compressed_shape))
            data = [num_samples] + shape
        else:
            data = [num_samples] + data
        shapeLibrary.insert(id, data)
    return shapeLibrary

def is_float(element: any) -> bool:
    try:
        float(element)
        return True
    except ValueError:
        return False
            
def findPeriodicPattern(listOfEvents, currentPattern):
    # If no pattern is found, returns 0
    # If the list is empty, returns -1

    if len(listOfEvents) == 0:
        return -1
    elif len(listOfEvents) == 1:
        return 0
    else:
        counter = 0
        for event in listOfEvents:
            if event == currentPattern:
                counter += 1
        period = int(len(listOfEvents) / counter)
        return period
    
def findVariableIndexes(variableEvents):
    variableIndexesList = []
    for event in variableEvents:
        for index in range(0, len(event)):
            if index not in variableIndexesList and event[index] != variableEvents[0][index]:
                variableIndexesList.append(index)
    return variableIndexesList

def evaluatePeriodicity(listOfEvents):
    periodicEvents = []
    variableEvents = []
    for event in listOfEvents:
        eventFound = False
        for periodicEvent in periodicEvents:
            if event == periodicEvent[0]:
                eventFound = True
        if eventFound == True:
            continue
        else:
            period = findPeriodicPattern(listOfEvents, event)
            if period != len(listOfEvents):
                periodicEvents.append([event, period])
            else:
                variableEvents.append(event)
                
    variableIndexes = findVariableIndexes(variableEvents)

    return periodicEvents, variableEvents, variableIndexes
    
def extractRFevents(seq, periodicEvents, variableEvents, variableIndexes):
    if 1 in variableIndexes:
        pass
    else:
        for event in periodicEvents:
            if event[0][1] != 0:
                rfFormat = seq.rfLibrary.data[event[0][1]]
                rfAmplitude = rfFormat[0]
                rfMagWaveform = seq.shapeLibrary.data[rfFormat[1]][1:]
                rfPhaseWaveform = seq.shapeLibrary.data[rfFormat[2]][1:]
                rfWaveform = []
                for i in range(0,len(rfMagWaveform)):
                    rfWaveform.append(rfMagWaveform[i])
                    rfWaveform.append(rfPhaseWaveform[i])
                rfTimeShape = seq.shapeLibrary.data[rfFormat[3]][1:]
                rfDelayUs = rfFormat[4]
                rfFrequencyOffset = rfFormat[5]
                rfPhaseOffset = rfFormat[6]

                # Calculating flip angle for SinC RF pulse using area under curve
                gamma = 42.577 # MHz/T gyromagnetic ratio
                rfMagWaveformAmpl = []
                for value in rfMagWaveform:
                    rfMagWaveformAmpl.append(value * rfAmplitude * 1.35) # TO DO check the 1.35 value, it should be 1.
                areaUnderCurve = simpson(rfMagWaveformAmpl, dx=seq.rfRasterTime)
                flipAngle = int(areaUnderCurve * gamma * 2 * np.pi) 
                print("flipAngle", flipAngle)

                # Calculating bandwith time product for SinC RF pulse
                rfMagWaveform = np.array(rfMagWaveform)
                # fftRfWaveform = np.fft.fftshift(np.abs(np.fft.fft(rfMagWaveform/rfMagWaveform.max()))) / np.sqrt(len(rfMagWaveform))
                # preciseFftRfWaveformEq = interpolate.InterpolatedUnivariateSpline(rfTimeShape, fftRfWaveform)
                # preciseFftRfWaveform = list(preciseFftRfWaveformEq(linspace(0, rfTimeShape[-1], 64000)))
                rfMagWaveform = list(rfMagWaveform)
                diffRfMagWaveform = np.diff(np.sign(np.diff(rfMagWaveform)))
                # plt.plot(rfMagWaveform)
                # plt.show()
                rfMinIndex = []
                for index in range(0, len(diffRfMagWaveform)):
                    if diffRfMagWaveform[index] > 0:
                        rfMinIndex.append(index)
                rfMaxIndex = []
                for index in range(0, len(diffRfMagWaveform)):
                    if diffRfMagWaveform[index] < 0:
                        rfMaxIndex.append(index)
                print("rfMinIndex", rfMinIndex)
                print("rfMaxIndex", rfMaxIndex)
                centralMaxIndex = int(len(rfMinIndex)/2)
                print("centralMaxIndex", centralMaxIndex)
                peakHalfWidth = (rfMaxIndex[centralMaxIndex] - rfMinIndex[centralMaxIndex-1]) * seq.rfRasterTime * 2
                print("peakHalfWidth", peakHalfWidth)
                lobeWidthTable =[]
                for index in range(0, len(rfMaxIndex)):
                    if index <= centralMaxIndex:
                        if index == 0:
                            lobeWidth = rfMinIndex[index] * seq.rfRasterTime * 2
                        else:
                            lobeWidth = (rfMinIndex[index] - rfMinIndex[index-1]) * seq.rfRasterTime * 2
                        lobeWidthTable.append(lobeWidth)
                print("lobeWidthTable", lobeWidthTable)
                peakHalfWidth= max(lobeWidthTable)/2
                # TO DO remove the *2 in bandWith and timeBwProduct when the original file is corrected
                bandWidth = 0.65 / (peakHalfWidth) 
                timeBwProduct = len(rfMagWaveform) * 2 * seq.rfRasterTime * bandWidth
                print("timeBwProduct", timeBwProduct)

                # Calculating slice thickness
                gamma = 42.577 # MHz/T gyromagnetic ratio
                sliceThickness = 5.0 # default value in mm
                if event[0][4] != 0:
                    sliceSelectionGradient = seq.gradLibrary.data[event[0][4]]
                    sliceSelectionGradientAmplitude = sliceSelectionGradient[0]
                    correctedAmplitude_mTm = (sliceSelectionGradientAmplitude * 1e-3) / gamma
                    sliceThickness = bandWidth / (gamma * correctedAmplitude_mTm) # in mm
                print("sliceThickness", sliceThickness)
                

                # sdlRfEvent = Rf()
                sdlRfExcitation = RfExcitation()
                sdlRfExcitation.array = rfWaveform
                # TO DO correct duration after the orignal file is corrected
                sdlRfExcitation.duration = len(rfMagWaveform) * 2 * seq.rfRasterTime * 1e6
                sdlRfExcitation.initial_phase = rfPhaseOffset
                sdlRfExcitation.freq_offset = rfFrequencyOffset
                sdlRfExcitation.thickness = sliceThickness
                # print(seq.timeBwProduct)

