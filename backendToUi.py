################################################################################
### mtrk project - Transition from backend to UI.                            ###
### Version 0.1.1                                                            ###
### Anais Artiges and the mtrk project team at NYU - 04/29/2024              ###
################################################################################  

import json
import jsbeautifier
import re
import ast
from pprint import pprint
from numpy import add
from sdlFileCreator import *
import os
import copy

#############################################################
### Creating SDL file from web-based UI
#############################################################
"""
Documentation for this module.

"""

def create_sdl_from_ui_inputs(block_to_box_objects, block_structure, 
                              block_to_loops, block_to_duration, 
                              block_number_to_block_object, configurations):
    """
    Create an SDL file from web-based UI inputs.

    Args:
        block_to_box_objects (dict): Mapping of block names to box objects.
        block_structure (dict): Mapping of block names to their structure.
        block_to_loops (dict): Mapping of block names to the number of loops.
        block_to_duration (dict): Mapping of block names to their duration.
        block_number_to_block_object (dict): Mapping of block numbers to block objects.
        configurations (dict): Configuration settings.

    Returns:
        None
    """
    ### Initialize SDL file
    ## TO DO - need to intialize without loading file
    file_path = os.path.abspath("mtrk_designer_api/init_data/miniflash.mtrk")
    with open(file_path) as sdlFile:
        sdlData = json.load(sdlFile)
        sequence_data = PulseSequence(**sdlData)
    sdlInitialize(sequence_data)

    sequence_data.file = File()
    sequence_data.infos = Info()
    sequence_data.settings = Settings()
    sequence_data.instructions = {}
    sequence_data.objects = {}
    sequence_data.arrays = {}
    sequence_data.equations = {}
    
    ### Special initialization for very simple sequences with no block structure
    if not block_to_loops or block_to_loops == {'Main': 1}:
        block_to_loops.update({"simple_main_block": 1})
        block_structure.update({"Main": ["simple_main_block"]})
        block_to_box_objects.update({"simple_main_block": copy.deepcopy(block_to_box_objects["Main"])})
        block_to_duration.update({"simple_main_block": 0})
        for index in range(len(block_to_box_objects["Main"])):
            block_to_box_objects["Main"][index]["type"] = "Block"
            block_to_box_objects["Main"][index]["start_time"] = 0
            block_to_box_objects["Main"][index]["block"] = 1
            block_to_box_objects["Main"][index]["name"] = "simple_main_block"

    updateSDLFile(sequence_data, block_to_box_objects, configurations,
                  block_number_to_block_object, block_to_loops, block_structure,
                  block_to_duration)
    
    ### writing of json schema to SDL file with formatting options
    with open('output.mtrk', 'w') as sdlFileOut:
        options = jsbeautifier.default_options()
        options.indent_size = 4
        data_to_print = jsbeautifier.beautify(\
                     json.dumps(sequence_data.model_dump(mode="json")), options)
        sdlFileOut.write(re.sub(r'}, {', '},\n            {', data_to_print)) 
        #purely aesthetic

def updateSDLFile(sequence_data, boxes, configurations, 
                  block_number_to_block_object, block_to_loops, block_structure,
                  block_to_duration):
    """
    Update the SDL file with new information.

    Args:
        sequence_data (PulseSequence): The SDL sequence data.
        boxes (dict): Mapping of box keys to box lists.
        configurations (dict): Configuration settings.
        block_number_to_block_object (dict): Mapping of block numbers to block objects.
        block_to_loops (dict): Mapping of block names to the number of loops.
        block_to_duration (dict): Mapping of block names to their duration.

    Returns:
        None
    """
    keys = boxes.keys()
    with open('test_backend.txt', 'w') as sdlFileOut:
        sdlFileOut.write("configurations \n") 
        sdlFileOut.write(str(configurations)) 
        sdlFileOut.write("\n\n") 
        sdlFileOut.write("block_number_to_block_object \n")
        sdlFileOut.write(str(block_number_to_block_object))
        sdlFileOut.write("\n\n") 
        sdlFileOut.write("block_to_loops \n") 
        sdlFileOut.write(str(block_to_loops))
        sdlFileOut.write("\n\n") 
        sdlFileOut.write("block_structure \n") 
        sdlFileOut.write(str(block_structure))
        sdlFileOut.write("\n\n") 
        sdlFileOut.write("block_to_duration \n") 
        sdlFileOut.write(str(block_to_duration)) 
        sdlFileOut.write("\n\n") 
        sdlFileOut.write("boxes \n") 
        sdlFileOut.write(str(boxes))
        sdlFileOut.write("\n\n") 
    
    instructionHeader = ["Default message", "off"]
    instructionEndTime = 0
    for boxKey in keys:
        boxList = boxes[boxKey]
        for box in boxList:
            if box["type"] == "event":
                box["type"] = box["axis"]
                box["axis"] = "event"
            if boxKey == "Main":
                instructionName = "main"
                instructionHeader = ["Running main loop", "off"]
            else:
                instructionName = boxKey
                if box["type"] != "Block":
                    specialTimingFlag = False
                    for value in block_number_to_block_object:
                        if block_number_to_block_object[value]["name"] == boxKey:
                            if block_number_to_block_object[value][\
                                                               "print_counter"] == True:
                                printCounter = "on"
                            else:
                                printCounter = "off"
                            instructionHeader = [block_number_to_block_object[value]["message"],
                                                 printCounter]
                        if block_number_to_block_object[value]["name"] == boxKey and \
                            block_number_to_block_object[value]["use_duration_equation"] == True:
                            equationEndTime = EquationRef()
                            equationEndTime.equation = block_number_to_block_object[value]["duration_equation_info"]["name"]
                            equation = block_number_to_block_object[value]["duration_equation_info"]["expression"]
                            instructionEndTime = [equation, equationEndTime]
                            specialTimingFlag = True
                        else:
                            if specialTimingFlag == False:
                                instructionEndTime = block_to_duration[boxKey]

            addInstruction(sequence_data, instructionName)
            instructionInformationList = getInstructionInformation(
                                            boxes = boxList,
                                            instructionName = instructionName,
                                            instructionHeader = instructionHeader,
                                            instructionEndTime = instructionEndTime)
            completeInstructionInformation(
                            sequence_data = sequence_data, 
                            instructionInformationList = instructionInformationList)

    makeBlockStructure(sequence_data, block_structure, block_to_loops)    
    
    fileInformationList = getFileInformation(configurations = configurations)
    completeFileInformation(sequence_data = sequence_data, 
                            fileInformationList = fileInformationList)
    settingsInformationList = getSequenceSettingsInformation(
                                                configurations = configurations)
    completeSequenceSettings(sequence_data = sequence_data, 
                             settingsInformationList = settingsInformationList)
    sequenceInfoInformationList = getSequenceInfoInformation(
                                                configurations = configurations)
    completeSequenceInformation(
                      sequence_data = sequence_data, 
                      sequenceInfoInformationList = sequenceInfoInformationList)

def makeBlockStructure(sequence_data, block_structure, block_to_loops):
    counterIndex = 1 
    for block in block_structure:
        if block == "Main":
            mainLoop = Loop(counter = counterIndex, range = 1, steps=[])   
            counterIndex += 1
            sequence_data.instructions["main"].steps.append(mainLoop)
            structureName = sequence_data.instructions["main"].steps[0]
        else:
            structureName = sequence_data.instructions[block]
        for blockName in block_structure[block]:
            if int(block_to_loops[blockName]) == 1:
                runBlock = RunBlock(block = blockName)
                structureName.steps.append(runBlock)
            else:
                loopBlock = Loop(counter = counterIndex, range = block_to_loops[blockName], steps=[])
                counterIndex += 1
                structureName.steps.append(loopBlock)
                runBlock = RunBlock(block = blockName)
                for step in structureName.steps:
                    if isinstance(step, Loop):
                        step.steps.append(runBlock)
        if block != "Main":
            submit_event = Submit()
            structureName.steps.append(submit_event)
            

#############################################################
### Functions to get new values from the web-based UI
#############################################################

def getFileInformation(configurations):
    """
    Get file information from the web-based UI.

    Args:
        configurations (dict): Configuration settings.

    Returns:
        list: List of file information.
    """
    formatInfo = configurations["file"]["format"]
    versionInfo = configurations["file"]["version"]
    measurementInfo = configurations["file"]["measurement"]
    systemInfo = configurations["file"]["system"]
    fileInformationList = [formatInfo, versionInfo, measurementInfo, 
                           systemInfo]
    return fileInformationList

def getSequenceSettingsInformation(configurations):
    """
    Get sequence settings information from the web-based UI.

    Args:
        configurations (dict): Configuration settings.

    Returns:
        list: List of sequence settings information.
    """
    readoutOsInfo = configurations["settings"]["readout"]
    variables = configurations["settings"]["variables"]
    settingsInformationList = [readoutOsInfo, variables]
    return settingsInformationList    

def getSequenceInfoInformation(configurations):
    """
    Get sequence info information from the web-based UI.

    Args:
        configurations (dict): Configuration settings.

    Returns:
        list: List of sequence info information.
    """
    descriptionInfo = configurations["info"]["description"]
    slicesInfo = configurations["info"]["slices"]
    fovInfo = configurations["info"]["fov"]
    # pelinesInfo = configurations["info"]["pelines"]
    pelinesInfo = configurations["info"]["resolution"]
    seqstringInfo = configurations["info"]["seqstring"]
    reconstructionInfo = configurations["info"]["reconstruction"]
    sequenceInfoInformationList = [descriptionInfo, slicesInfo, fovInfo, 
                                   pelinesInfo, seqstringInfo, 
                                   reconstructionInfo]
    return sequenceInfoInformationList

def getInstructionInformation(boxes, instructionName, instructionHeader, instructionEndTime):
    """
    Get instruction information from the web-based UI.

    Args:
        boxes (list): List of box objects.
        instructionName (str): Name of the instruction.
        instructionHeader (list): List containing the instruction header information.

    Returns:
        list: List of instruction information.
    """
    printMessageInfo = instructionHeader[0]
    printCounterInfo = instructionHeader[1]
    allStepInformationLists = []
    initInfo = ""
    for box in boxes:
        stepInformationList = getStepInformation(box)
        if box["type"] == "loop" and stepInformationList in \
                                                        allStepInformationLists:
            pass
        elif box["type"] == "init":
            initInfo = box["inputInitActionGradients"]
        else:
            allStepInformationLists.append(stepInformationList)
    instructionInformationList = [instructionName, printMessageInfo,
                                  printCounterInfo, instructionEndTime,
                                  allStepInformationLists, initInfo]
    return instructionInformationList
            
def getStepInformation(box):
    """
    Get step information from the web-based UI.

    Args:
        box (dict): Box object.

    Returns:
        list: List of step information.
    """
    actionName = box["type"]
    stepInformationList = [actionName]
    match actionName:
        case "run_block":
            blockName = box["name"]
            stepInformationList.extend([blockName])

        case "loop":
            counterInfo = box["block"]
            rangeInfo = box["loop_number"]
            ## TO DO decide either giving directly the all_step_info_lists or
            ## giving a list of step information from which it is extracted.
            derivedBox = {'type': 'run_block', 'name': box["name"]} 
            allStepInformationLists = []
            # for stepInformationList in allStepInformationLists:
            newStepInformationList = getStepInformation(derivedBox)
            allStepInformationLists.append(newStepInformationList)
            stepInformationList.extend([counterInfo, rangeInfo, 
                                        allStepInformationLists])  
            
        case "calc":
            typeInfo = box["inputCalcActionType"]
            floatInfo = box["inputCalcFloat"]
            incrementInfo = box["inputCalcIncrement"]
            stepInformationList.extend([typeInfo, floatInfo, incrementInfo]) 

        case "init":
            pass

        case "sync":
            objectInfo = box["inputSyncObject"]
            objectInformationList = getObjectInformation(typeInfo = actionName, 
                                                         box = box)
            timeInfo = int(float(box["inputSyncTime"]))
            stepInformationList.extend([objectInfo, objectInformationList,
                                        timeInfo]) 
            
        case "grad":
            axisInfo = box["axis"]
            objectInfo = box["name"]
            objectInformationList = getObjectInformation(typeInfo = actionName, 
                                                         box = box)
            stepInformationList.extend([axisInfo, objectInfo, 
                                        objectInformationList])
            
            if box["use_equation_time"] == True:
                timeTypeInfo = "equation"
                timeEquationNameInfo = box["equation_time_info"]["name"]
                stepInformationList.extend([timeTypeInfo, 
                                            timeEquationNameInfo])
                equationTimeInfo = box["equation_time_info"]["expression"]
                stepInformationList.extend([equationTimeInfo])
            else:
                timeInfo = int(float(box["start_time"]))
                stepInformationList.extend([timeInfo])
 
            if str(box["flip_amplitude"]) == "True":
                flipAmplitudeInfo = "flip"
                stepInformationList.extend([flipAmplitudeInfo])
            if str(box["variable_amplitude"]) == "True":
                amplitudeTypeInfo = "equation"
                amplitudeEquationNameInfo = box["equation_info"]["name"]
                stepInformationList.extend([amplitudeTypeInfo, 
                                            amplitudeEquationNameInfo])
                equationInfo = box["equation_info"]["expression"]
                stepInformationList.extend([equationInfo])

        case "rf":
            objectInfo= box["name"]
            objectInformationList = getObjectInformation(typeInfo = actionName, 
                                                         box = box)
            stepInformationList.extend([objectInfo, objectInformationList])
            if box["use_equation_time"] == True:
                timeTypeInfo = "equation"
                timeEquationNameInfo = box["equation_time_info"]["name"]
                stepInformationList.extend([timeTypeInfo, 
                                            timeEquationNameInfo])
                equationTimeInfo = box["equation_time_info"]["expression"]
                stepInformationList.extend([equationTimeInfo])
            else:
                timeInfo = int(float(box["start_time"]))
                stepInformationList.extend([timeInfo])
                
            if box["use_equation_freqoffset"] == True:
                timeTypeInfo = "equation"
                timeEquationNameInfo = box["equation_freqoffset_info"]["name"]
                stepInformationList.extend([timeTypeInfo, 
                                            timeEquationNameInfo])
                equationFreqOffsetInfo = box["equation_freqoffset_info"]["expression"]
                stepInformationList.extend([equationFreqOffsetInfo])
            else:
                freqoffsetInfo = int(float(box["freq_offset"]))
                stepInformationList.extend([freqoffsetInfo])
            addedPhaseTypeInfo = box["rf_added_phase_type"]
            addedPhaseFloatInfo = box["rf_added_phase_float"]
            stepInformationList.extend([addedPhaseTypeInfo, 
                                        addedPhaseFloatInfo])
        
        case "adc":
            objectInfo = box["name"]
            objectInformationList = getObjectInformation(typeInfo = actionName, 
                                                         box = box)
            stepInformationList.extend([objectInfo, objectInformationList])

            if box["use_equation_time"] == True:
                timeTypeInfo = "equation"
                timeEquationNameInfo = box["equation_time_info"]["name"]
                stepInformationList.extend([timeTypeInfo, 
                                            timeEquationNameInfo])
                equationTimeInfo = box["equation_time_info"]["expression"]
                stepInformationList.extend([equationTimeInfo])
            else:
                timeInfo = int(float(box["start_time"]))
                stepInformationList.extend([timeInfo])

            frequencyInfo = box["frequency"]
            phaseInfo = box["phase"]
            addedPhaseTypeInfo = box["adc_added_phase_type"]
            addedPhaseFloatInfo = box["adc_added_phase_float"]
            mdhInfoList = json.loads(box["mdh"])
            # print("mdhInfoList ", mdhInfoList)
            stepInformationList.extend([frequencyInfo, phaseInfo,
                                        addedPhaseTypeInfo,addedPhaseFloatInfo,
                                        mdhInfoList])
            
        # case "mark":
        #     print("+-+-+ mark box[equation_time_info]", box["equation_time_info"])
        #     if box["use_equation_time"] == True:
        #         timeTypeInfo = "equation"
        #         timeEquationNameInfo = box["equation_time_info"]["name"]
        #         stepInformationList.extend([timeTypeInfo, 
        #                                     timeEquationNameInfo])
        #         equationTimeInfo = box["equation_time_info"]["expression"]
        #         stepInformationList.extend([equationTimeInfo])
        #     else:
        #         timeInfo = int(float(box["start_time"]))
        #         stepInformationList.extend([timeInfo])

        case "submit":
            pass
        
        case _:
            # print("The type " + actionName + " could not be identified.")
            pass
    return stepInformationList
            
def getObjectInformation(typeInfo, box):
    """
    Get the information of an object based on its type.

    Args:
        typeInfo (str): The type of the object.
        box (dict): The box containing the object information.

    Returns:
        list: A list containing the object information.

    Raises:
        None

    """
    objectInformationList = [typeInfo]
    match typeInfo:
        case "rf":
            ## TO DO make the duration step flexible
            durationInfo = len(box["array_info"]["array"])*20
            arrayInfo = box["array_info"]["name"]
            arrayInformationList = getArrayInformation(box = box)
            initPhaseInfo = box["init_phase"]
            freqOffsetInfo = box["freq_offset"]
            thicknessInfo = box["thickness"]
            flipAngleInfo = box["flip_angle"]
            purposeInfo = box["purpose"]
            objectInformationList.extend([durationInfo, arrayInfo, 
                                          arrayInformationList, 
                                          initPhaseInfo, freqOffsetInfo, thicknessInfo, 
                                          flipAngleInfo, purposeInfo])
            
        case "grad":
            durationInfo = len(box["array_info"]["array"])*10
            arrayName = box["array_info"]["name"]
            arrayInformationList = getArrayInformation(box = box)
            # TO DO add "tail" to the dictionnary
            # tailInfo = box["tail"]
            tailInfo = 0
            amplitudeInfo = box['amplitude']
            objectInformationList.extend([durationInfo, arrayName, 
                                          arrayInformationList, 
                                          tailInfo, amplitudeInfo])
            
        case "adc":
            durationInfo = box["adc_duration"]*1e3
            samplesInfo = box["samples"]
            dwelltimeInfo = int(box["dwell_time"])
            objectInformationList.extend([durationInfo, samplesInfo, 
                                          dwelltimeInfo])

        case "sync":
            eventInfo = box["inputSyncEventParam"]
            durationInfo = box["inputSyncDuration"]
            objectInformationList.extend([durationInfo, eventInfo])

        case "init":
            pass

        case "calc":
            pass

        case _:
            # print("The type " + typeInfo + " could not be identified.")
            pass

    return objectInformationList
            
def getArrayInformation(box):
    """
    Retrieves information about an array from a given box.

    Args:
        box (dict): The box containing the array information.

    Returns:
        list: A list containing the encoding, type, size, and data of the array.
    """
    # TO DO add "encoding" to the dictionnary
    # encodingInfo = box["encoding"]
    encodingInfo = "text"
    # TO DO add "type" to the dictionnary
    # typeInfo = box["type"]
    if box["type"] == "rf":
        typeInfo = "complex_float"
        magn_array = box["array_info"]["array"]
        phase_array = box["phase_array_info"]["array"]
        array = []
        if len(magn_array) == len(phase_array):
            for index in range(0,len(magn_array)):
                array.append(magn_array[index])
                array.append(phase_array[index]) 
        sizeInfo = len(array)/2
    else:
        typeInfo = "float"
        array = box["array_info"]["array"]
        sizeInfo = len(array)
    dataInfoList = array

    arrayInformationList = [encodingInfo, typeInfo, sizeInfo, dataInfoList]
    return arrayInformationList
