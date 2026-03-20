################################################################################
### mtrk project - Console interface to create an SDL file.                  ###
### Version 0.1.1                                                            ###
### Anais Artiges and the mtrk project team at NYU - 09/07/2023              ###
################################################################################  

from SDL_read_write.pydanticSDLHandler import *
from sdlFileCreator import *

#############################################################
### User interface interaction
#############################################################

def mtrkConsoleUI(sequence_data):
    """
    This function provides a console-based user interface for configuring a sequence.

    Args:
        sequence_data (SequenceData): The sequence data object to be modified.

    Returns:
        SequenceData: The modified sequence data object.
    """
    sdlInitialize(sequence_data)

    ### file section
    print("*************** - FILE - ***************")
    print("Do you want to provide file information? (yes/no)")
    if(input() == "yes"):
        fileInformationList = getFileInformation()
        completeFileInformation(sequence_data, fileInformationList)
    else:
        print("Default File information used.")

    ### settings section
    print("*************** - SETTINGS - ***************")
    print("Do you want to provide general sequence settings? (yes/no)")
    if(input() == "yes"):
        settingsInformationList = getSequenceSettings()
        completeSequenceSettings(sequence_data, settingsInformationList)
    else:
        print("Default Settings information used.")

    ### info section
    print("*************** - INFORMATION - ***************")
    print("Do you want to provide general sequence information? (yes/no)")
    if(input() == "yes"):
        sequenceInfoInformationList = getSequenceInformation()
        completeSequenceInformation(sequence_data, sequenceInfoInformationList)
    else:
        print("Default Info information used.")

    ### instructions section
    print("*************** - INSTRUCTIONS - ***************")
    answer = "yes"
    while(answer == "yes"):
        print("*** Do you want to add a new instruction? (yes/no)")
        answer = input()
        if(answer == "yes"):
            print("Instruction name (str):")
            instructionName = input()
            addInstruction(sequence_data, instructionName)
            print("Do you want to provide instruction information? (yes/no)")
            if(input() == "yes"):
                instructionInformationList = \
                                      getInstructionInformation(instructionName)
                completeInstructionInformation(sequence_data, 
                                               instructionInformationList)
                print("instructionToModify = "+ \
                      str(sequence_data.instructions[instructionName]))
            else:
                print("Default Instruction information used.")
        else:
            print("*******************************************")

    return sequence_data

#############################################################
### Functions to get new values from the user
#############################################################

def getFileInformation():
    """
    Prompts the user to enter file information and returns a list containing the entered information.

    Returns:
        list: A list containing the format, version, measurement, and system information entered by the user.
    """
    print("format (str)")
    formatInfo = input()
    print("version (int)")
    versionInfo = int(input())
    print("measurement (str)")
    measurementInfo = input()
    print("system (str)")
    systemInfo = input()
    fileInformationList = [formatInfo, versionInfo, measurementInfo, systemInfo]
    return fileInformationList

def getSequenceSettings():
    """
    Retrieves the sequence settings from the user.

    Returns:
        list: A list containing the sequence settings information.
    """
    print("readout_os (int)")
    readoutOsInfo = input()
    settingsInformationList = [readoutOsInfo]
    return settingsInformationList

def getSequenceInformation():
    """
    Prompts the user to enter sequence information and returns a list of the entered information.

    Returns:
        list: A list containing the following sequence information:
            - description (str)
            - slices (int)
            - fov (int)
            - pelines (int)
            - seqstring (str)
            - reconstruction (str)
    """
    print("description (str)")
    descriptionInfo = input()
    print("slices (int)")
    slicesInfo = int(input())
    print("fov (int)")
    fovInfo = int(input())
    print("pelines (int)")
    pelinesInfo = int(input())
    print("seqstring (str)")
    seqstringInfo = input()
    print("reconstruction (str)")
    reconstructionInfo = input()
    sequenceInfoInformationList = [descriptionInfo, slicesInfo, fovInfo, 
                                   pelinesInfo, seqstringInfo, 
                                   reconstructionInfo]
    return sequenceInfoInformationList

def getInstructionInformation(instructionName):
    """
    Retrieves information for a specific instruction.

    Args:
        instructionName (str): The name of the instruction.

    Returns:
        list: A list containing the instruction name, print message info, print counter info, and a list of step information lists.
    """
    print("Message to print (str): ")
    printMessageInfo = input()
    print("Printing counter option (on/off): ")
    printCounterInfo = input()
    stepAnswer = "yes"
    stepIndex = 0
    allStepInformationLists = []
    while(stepAnswer == "yes"):
        print("Do you want to add a new step? (yes/no)")
        stepAnswer = input()
        if(stepAnswer == "yes"):
            print("Do you want to provide step information? (yes/no)")
            if(input() == "yes"):
                newStepInformationList = getStepInformation()
                allStepInformationLists.append(newStepInformationList)
            else:
                print("Default Instruction information used.")
            stepIndex += 1
        else:
            pass
    instructionInformationList = [instructionName, printMessageInfo,
                                  printCounterInfo, allStepInformationLists]
    return instructionInformationList

def getStepInformation():
    """
    Prompts the user to provide step information based on the selected action type.
    Returns a list containing the step information.

    Returns:
        list: A list containing the step information based on the selected action type.
    """
    print("Provide step action type: ")
    print("Action (run_block/loop/calc/init/sync/grad/rf/adc/mark/submit): ")
    actionName = input()
    print("Provide information for step of type " + str(actionName) + ": ")
    stepInformationList = [actionName]
    match actionName:
        case "run_block":
            print("block (str)")
            blockName = input()
            print("Do you want to provide block information? (yes/no)")
            blockInformationList = []
            if(input()=="yes"):
                blockInformationList = getInstructionInformation(blockName)
            else:
                print("Default Array information used.")
            stepInformationList.extend([blockName, blockInformationList])
        case "loop":
            print("counter (int)")
            counterInfo = input()
            print("range (int)")
            rangeInfo = input()
            print("steps (Step)")
            stepAnswerLoop = "yes"
            stepIndexLoop = 0
            allStepInformationLists = []
            while(stepAnswerLoop == "yes"):
                print("Do you want to add a new step in the loop? (yes/no)")
                stepAnswerLoop = input()
                nextAnswer = "yes"
                if(stepAnswerLoop == "yes"):
                    print("Do you want to provide step information? (yes/no)")
                    nextAnswer = input()
                    if(nextAnswer == "yes"):
                        newStepInformationList = getStepInformation()
                        allStepInformationLists.append(newStepInformationList)
                    else:
                        print("Default Loop information used.")
                    stepIndexLoop += 1
                else:
                    pass  
            stepInformationList.extend([counterInfo, rangeInfo, 
                                        allStepInformationLists])  
        case "calc":
            print("type (str)")
            typeInfo = input()
            print("float (float)")
            floatInfo = input()
            print("increment (int)")
            incrementInfo = input() 
            stepInformationList.extend([typeInfo, floatInfo, incrementInfo]) 
        case "init":
            print("gradients (str)")
            gradientInfo = input()
            stepInformationList.extend([gradientInfo]) 
        case "sync":
            print("object (str)")
            objectInfo = input()
            print("Do you want to provide object information? (yes/no)")
            if(input()=="yes"):
                objectInformationList = getObjectInformation(actionName)
            else:
                print("Default Object information used.")
            print("time (int)")
            timeInfo = input()
            stepInformationList.extend([objectInfo, objectInformationList, 
                                        timeInfo]) 
        case "grad":
            print("axis (slice/read/phase)")
            axisInfo = input()
            print("object (str)")
            objectInfo = input()
            print("Do you want to provide object information? (yes/no)")
            if(input()=="yes"):
                objectInformationList = getObjectInformation(actionName)
            else:
                print("Default Object information used.")
            print("time (int)")
            timeInfo = input() 
            stepInformationList.extend([axisInfo, objectInfo, 
                                        objectInformationList, timeInfo]) 
            print("Do you want to add an amplitude option? (yes/no)")
            if(input()=="yes"):
                print("amplitude option (flip/equation)")
                amplitudeAnswer = input()
                if(amplitudeAnswer=="flip"):
                    flipAmplitudeInfo = "flip"
                    stepInformationList.extend([flipAmplitudeInfo])
                elif(amplitudeAnswer=="equation"):
                    amplitudeTypeInfo = "equation"
                    print("amplitude equation name (str)")
                    amplitudeEquationNameInfo = input()
                    stepInformationList.extend([amplitudeTypeInfo, 
                                                amplitudeEquationNameInfo])
                    print("Do you want to complete equation information? (yes/no)")
                    if(input()=="yes"):
                        print("equation (str)")
                        equationInfo = input()
                        stepInformationList.extend([equationInfo])
            else:
                print("No amplitude option added.")

        case "rf":
            print("object (str)")
            objectInfo = input()
            print("Do you want to provide object information? (yes/no)")
            if(input()=="yes"):
                objectInformationList = getObjectInformation(actionName)
            else:
                print("Default Object information used.")
            print("time (float)")
            timeInfo = input()
            print("added_phase type (str)")
            addedPhaseTypeInfo = input()
            print("added_phase float (float)")
            addedPhaseFloatInfo = input()
            stepInformationList.extend([objectInfo, objectInformationList, 
                                        timeInfo, addedPhaseTypeInfo, 
                                        addedPhaseFloatInfo])
        case "adc":
            print("object (str)")
            objectInfo = input()
            print("Do you want to provide object information? (yes/no)")
            if(input()=="yes"):
                objectInformationList = getObjectInformation(actionName)
            else:
                print("Default Object information used.")
            print("time (float)")
            timeInfo = input()
            print("frequency (int)")
            frequencyInfo = input()
            print("phase (int)")
            phaseInfo = input()
            print("added_phase type (str)")
            addedPhaseTypeInfo = input()
            print("added_phase float (float)")
            addedPhaseFloatInfo = input()
            print("mdh (dict[str, MdhOption])")
            mdhInfoList = []
            #### TO DO !!! complete mdhInfoList
            print("passed for now") 
            stepInformationList.extend([objectInfo, objectInformationList, 
                                        timeInfo, frequencyInfo, phaseInfo,
                                        addedPhaseTypeInfo,addedPhaseFloatInfo,
                                        mdhInfoList])
        case "mark":
            print("time (float)")
            timeInfo = input()
            stepInformationList.extend([timeInfo])
        case "submit":
            print("Nothing to customize in submit.")
        case _:
            print("The type " + actionName + " could not be identified.")
    return stepInformationList

def getObjectInformation(typeInfo):
    """
    Get information for an object of the given type.

    Args:
        typeInfo (str): The type of the object.

    Returns:
        list: A list containing the object information based on the type.

    Raises:
        None

    """
    print("Provide information for object of type " + str(typeInfo) + ": ")
    print("duration (int)")
    durationInfo = input()
    objectInformationList = [typeInfo, durationInfo]
    match typeInfo:
        case "rf":
            print("array (str)")
            arrayInfo = input()
            print("Do you want to provide array information? (yes/no)")
            if(input()=="yes"):
                arrayInformationList = getArrayInformation()
            else:
                print("Default Array information used.")
            print("initial_phase (int)")
            initPhaseInfo = input()
            print("freq_offset (int)")
            freqOffsetInfo = input()
            print("thickness (int)")
            thicknessInfo = input()
            print("flipangle (int)")
            flipAngleInfo = input()
            print("purpose (str)")
            purposeInfo = input()
            objectInformationList.extend([arrayInfo, arrayInformationList, 
                                          initPhaseInfo, freqOffsetInfo, thicknessInfo, 
                                          flipAngleInfo, purposeInfo])
        case "grad":
            print("array (str)")
            arrayName = input()
            print("Do you want to provide array information? (yes/no)")
            if(input()=="yes"):
                arrayInformationList = getArrayInformation()
            else:
                print("Default Array information used.")
            print("tail (int)")
            tailInfo = input()
            print("amplitude (float)")
            amplitudeInfo = input()
            objectInformationList.extend([arrayName, arrayInformationList, 
                                          tailInfo, amplitudeInfo])
        case "adc":
            print("samples (int)")
            samplesInfo = input()
            print("dwelltime (int)")
            dwelltimeInfo = input()
            objectInformationList.extend([samplesInfo, dwelltimeInfo])
        case "sync":
            print("event (str)")
            eventInfo = input()
            objectInformationList.extend([eventInfo])
        case _:
            print("The type " + typeInfo + " could not be identified.")
    return objectInformationList

def getArrayInformation():
    """
    Prompts the user to enter information about an array and returns the array information as a list.

    Returns:
        list: A list containing the following information about the array:
            - encoding (str)
            - type (str)
            - size (int)
            - data (list of floats)
    """
    print("encoding (str)")
    encodingInfo = input()
    print("type (str)")
    typeInfo = input()
    print("size (int)")
    sizeInfo = input()
    print("data (float, float, ...)")
    dataInfoList = [float(elem) for elem in input().split(", ")]
    arrayInformationList = [encodingInfo, typeInfo, sizeInfo, dataInfoList]
    return arrayInformationList
