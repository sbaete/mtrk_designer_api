################################################################################
### mtrk project - SDL file generator allowing to simply define MRI pulse    ###
### sequences that can be read by the mtrk project simulator and driver      ###
### sequence.                                                                ###
### Version 0.1.1                                                            ###
### Anais Artiges and the mtrk project team at NYU - 04/29/2024              ###
################################################################################  

from SDL_read_write.pydanticSDLHandler import *

#############################################################
### Functions to create SDL file objects
#############################################################

## SDL file initialization with mandatory sections and default values
def sdlInitialize(sequence_data):
    """
    Initializes the sequence_data object with the necessary attributes.

    Args:
        sequence_data: The sequence_data object to be initialized.

    Returns:
        None
    """
    sequence_data.file = File()
    sequence_data.infos = Info()
    sequence_data.settings = Settings()
    sequence_data.instructions = {}
    sequence_data.objects = {}
    sequence_data.arrays = {}
    sequence_data.equations = {}


def addInstruction(sequence_data, instructionName):
    """
    Adds an instruction to the sequence data.

    Args:
        sequence_data (SequenceData): The sequence data object.
        instructionName (str): The name of the instruction.

    Returns:
        None
    """
    sequence_data.instructions[instructionName] = Instruction(steps=[])


def addStep(instructionToModify, stepIndex, actionName):
    """
    Adds a step to the given instruction.

    Args:
        instructionToModify (Instruction): The instruction object to modify.
        stepIndex (int): The index of the step to add.
        actionName (str): The name of the action for the step.

    Raises:
        ValueError: If the actionName is not available.

    Returns:
        None
    """
    try:
        match actionName:
            case "run_block":
                instructionToModify.steps.append(RunBlock())
                print("Disclaimer: no actual block created here, it will be \
                      created if you provide information for the step.")
            case "loop":
                instructionToModify.steps.append(Loop(steps=[]))
            case "calc":
                instructionToModify.steps.append(Calc())
            case "sync":
                instructionToModify.steps.append(Sync())
            case "grad":
                instructionToModify.steps.append(Grad())
            case "rf":
                instructionToModify.steps.append( 
                                                 Rf(added_phase = AddedPhase()))
            case "adc":
                instructionToModify.steps.append( 
                                        Adc(added_phase = AddedPhase(), mdh={}))
                # mdhOptAnswer = "yes"
                # while(mdhOptAnswer == "yes"):
                #     print("Do you want to add a new MDH option? (yes/no)")
                #     mdhOptAnswer = input()
                #     if(mdhOptAnswer == "yes"):
                #         addMdhOption(instructionToModify.steps, stepIndex)
                #     else:
                #         pass
            case "mark":
                instructionToModify.steps.append(Mark())
            case "submit":
                instructionToModify.steps.append(Submit())
            case _:
                raise ValueError(actionName + " is not available")
    except ValueError as e:
        print(str(e))


def addMdhOption(stepToModify, stepIndex):
    """
    Adds an MDH option to the specified step.

    Parameters:
    - stepToModify (list): The list of steps to modify.
    - stepIndex (int): The index of the step to modify.

    Returns:
    None
    """
    print("Provide MDH option information: ")
    print("MDH option type (str): ")
    stepToModify[stepIndex].mdh[input()] = MdhOption()


def addObject(sequence_data, objectName, typeName):
    """
    Add an object to the sequence_data based on the provided typeName.

    Args:
        sequence_data (SequenceData): The sequence data object to add the object to.
        objectName (str): The name of the object to add.
        typeName (str): The type of the object to add.

    Raises:
        ValueError: If the typeName is not available.

    Returns:
        None
    """
    try:
        match typeName:
            case "rf":
                sequence_data.objects[objectName] = RfExcitation()
            case "grad":
                sequence_data.objects[objectName] = GradientObject()
            case "adc":
                sequence_data.objects[objectName] = AdcReadout()
            case "sync":
                sequence_data.objects[objectName] = Ttl()
            case _:
                raise ValueError(typeName + " is not available")
    except ValueError as e:
        print(str(e))


def addArray(sequence_data, arrayName):
    """
    Adds an array to the sequence data.

    Parameters:
    - sequence_data: The sequence data object to add the array to.
    - arrayName: The name of the array to be added.

    Returns:
    None
    """
    sequence_data.arrays[arrayName] = GradientTemplate()


def addEquation(sequence_data, equationName):
    """
    Adds an equation to the sequence data.

    Parameters:
    - sequence_data (SequenceData): The sequence data object.
    - equationName (str): The name of the equation.

    Returns:
    None
    """
    sequence_data.equations[equationName] = Equation()


#############################################################
### Functions to fill SDL file objects with new values
#############################################################

def completeFileInformation(sequence_data, fileInformationList):
    """
    Completes the file information in the sequence_data object using the provided fileInformationList.

    Args:
        sequence_data (SequenceData): The sequence_data object to update.
        fileInformationList (list): A list containing the file information in the following order:
                                    [formatInfo, versionInfo, measurementInfo, systemInfo]

    Returns:
        None
    """
    ## fileInformationList = [formatInfo, versionInfo, measurementInfo, 
    ##                        systemInfo]
    if(fileInformationList != []):
        sequence_data.file = File(format = fileInformationList[0], 
                                  version = int(fileInformationList[1]), 
                                  measurement = fileInformationList[2], 
                                  system = fileInformationList[3])

def completeSequenceSettings(sequence_data, settingsInformationList):
    """
    Complete the sequence settings based on the provided information.

    Args:
        sequence_data (SequenceData): The sequence data object to update.
        settingsInformationList (list): A list of settings information.

    Returns:
        None
    """
    ## settingsInformationList = [readoutOsInfo]
    if(settingsInformationList != []):
        sequence_data.settings = Settings(readout_os = \
                                                     settingsInformationList[0])
    for variable in settingsInformationList[1]:
        if variable not in sequence_data.settings.model_fields:
            sequence_data.settings.set(variable, settingsInformationList[1][variable])
        else:
            sequence_data.settings.set(variable, settingsInformationList[1][variable])

def completeSequenceInformation(sequence_data, sequenceInfoInformationList):
    """
    Completes the sequence settings based on the provided information.

    Args:
        sequence_data (SequenceData): The sequence data object to update.
        settingsInformationList (list): A list of settings information.

    Returns:
        None
    """
    ## sequenceInfoInformationList = [descriptionInfo, slicesInfo, fovInfo, 
    ##                                pelinesInfo, seqstringInfo, 
    ##                                reconstructionInfo]
    if(sequenceInfoInformationList != []):
        sequence_data.infos = Info(
                               description = sequenceInfoInformationList[0], 
                               slices = int(sequenceInfoInformationList[1]), 
                               fov = int(sequenceInfoInformationList[2]),  
                               pelines = int(sequenceInfoInformationList[3]), 
                               seqstring = sequenceInfoInformationList[4],
                               reconstruction = sequenceInfoInformationList[5])

def completeInstructionInformation(sequence_data, instructionInformationList):
    """
    Updates the instruction information based on the provided list.

    Args:
        sequence_data (SequenceData): The sequence data object.
        instructionInformationList (list): A list containing the instruction information.
            The list should have the following structure:
            [instructionName, printMessageInfo, printCounterInfo, instructionEndTime, allStepInformationLists]

    Returns:
        None
    """
    ## instructionInformationList = [instructionName, printMessageInfo, 
    ##                               printCounterInfo, instructionEndTime,
    ##                               allStepInformationLists, initInfo]
    if(instructionInformationList != []):
        instructionToModify = \
                       sequence_data.instructions[instructionInformationList[0]]
        instructionToModify.print_message = instructionInformationList[1]
        printCounterOption = instructionInformationList[2]
        if(printCounterOption=="on" or printCounterOption=="off"):
            instructionToModify.print_counter = printCounterOption
        else:
            print(printCounterOption + 
                  " is not valid. It should be 'on' or 'off'.")
        allStepInformationLists = []
        for instruction in instructionInformationList[4]:
            if instruction != ['Block']:
                allStepInformationLists.append(instruction)
        
        if instructionInformationList[0] != "main":
            init_event = Init(gradients = instructionInformationList[5])
            instructionToModify.steps.append(init_event)


        for stepIndex in range(0, len(allStepInformationLists)):
            ## stepInformationList = [actionName, actionSpecificElements...]
            addStep(instructionToModify = instructionToModify, 
                    stepIndex = stepIndex + 1, 
                    actionName = allStepInformationLists[stepIndex][0])
            if instructionToModify.steps[0].action != "init":
                stepToModify = instructionToModify.steps[stepIndex]
            else:
                stepToModify = instructionToModify.steps[stepIndex + 1]
            completeStepInformation(sequence_data = sequence_data, 
                                    stepToModify = stepToModify, 
                                    stepInformationList = \
                                       allStepInformationLists[stepIndex])
        if instructionInformationList[0] != "main" and len(instructionToModify.steps) != 1:
            if type(instructionInformationList[3])!=int and type(instructionInformationList[3])!=float:
                mark_event = Mark(time = instructionInformationList[3][1])
                addEquation(sequence_data = sequence_data, 
                            equationName = mark_event.time.equation)
                sequence_data.equations[\
                        mark_event.time.equation].equation = \
                                                  instructionInformationList[3][0]
            else:
                mark_event = Mark(time = int(instructionInformationList[3]*100)*10)
            instructionToModify.steps.append(mark_event)
            submit_event = Submit()
            instructionToModify.steps.append(submit_event)

def completeStepInformation(sequence_data, stepToModify, stepInformationList):
    """
    Completes the step information based on the provided step information list.

    Args:
        sequence_data (SequenceData): The sequence data object.
        stepToModify (Step): The step object to modify.
        stepInformationList (list): The list containing the step information.

    Returns:
        None
    """
    ## stepInformationList = [actionName, actionSpecificElements...]
    if(stepInformationList != []):
        actionInfo = stepInformationList[0]
        match actionInfo:
            case "run_block":
                stepToModify.block = stepInformationList[1]
                ## stepInformationList = [actionName, blockName, 
                ##                        blockInformationList]
                ## NEEDED for console UI
                # stepToModify.block= stepInformationList[1]
                # addInstruction(sequence_data, stepToModify.block)
                # if(stepInformationList[2]!= []):
                #     completeInstructionInformation(
                #                                sequence_data = sequence_data, 
                #                                instructionInformationList = \
                #                                        stepInformationList[2])
                # else:
                #     print("Default Array information used.")
            case "loop":
                ## stepInformationList = [actionName, counterInfo, rangeInfo, 
                ##                        allStepInformationLists]
                savedStepToModify = stepToModify
                savedStepToModify.counter = stepInformationList[1]
                savedStepToModify.range= stepInformationList[2]
                for stepIndexLoop in range(0, len(stepInformationList[3])):
                    addStep(instructionToModify = savedStepToModify, 
                            stepIndex = stepIndexLoop, 
                            actionName = \
                                       stepInformationList[3][stepIndexLoop][0])
                    completeStepInformation(sequence_data = sequence_data, 
                                            stepToModify = \
                                         savedStepToModify.steps[stepIndexLoop], 
                                            stepInformationList = \
                                          stepInformationList[3][stepIndexLoop]) 
                stepToModify = savedStepToModify  
            case "calc":
                ## stepInformationList = [actionName, typeInfo, floatInfo, 
                ##                        incrementInfo]
                stepToModify.type = stepInformationList[1]
                stepToModify.float= stepInformationList[2]
                stepToModify.increment = stepInformationList[3] 
            # case "init":
            #     ## stepInformationList = [actionName, gradientInfo]
            #     stepToModify.gradients = stepInformationList[1]
            case "sync":
                ## stepInformationList = [actionName, objectInfo, 
                ##                        objectInformationList, timeInfo]
                stepToModify.object = stepInformationList[1]
                if stepToModify.object not in sequence_data.objects:
                    addObject(sequence_data=sequence_data, 
                            objectName=stepToModify.object,
                            typeName = "sync")
                completeObjectInformation(sequence_data = sequence_data, 
                                          objectName = stepToModify.object,
                                          objectInformationList = \
                                                         stepInformationList[2])
                stepToModify.time = stepInformationList[3]
            case "grad":
                ## stepInformationList = [actionName, axisInfo, objectInfo, 
                ##                        objectInformationList, timeInfo]
                stepToModify.axis = stepInformationList[1]
                stepToModify.object= stepInformationList[2]
                if stepToModify.object not in sequence_data.objects:
                    addObject(sequence_data=sequence_data, 
                              objectName=stepToModify.object,
                              typeName = "grad")
                completeObjectInformation(sequence_data = sequence_data, 
                                          objectName = stepToModify.object,
                                          objectInformationList = \
                                                         stepInformationList[3])
                
                if stepInformationList[4] == "equation":
                    stepToModify.time = EquationRef()
                    stepToModify.time.type = stepInformationList[4]
                    stepToModify.time.equation = stepInformationList[5]
                    addEquation(sequence_data = sequence_data, 
                            equationName = stepToModify.time.equation)
                    sequence_data.equations[\
                            stepToModify.time.equation].equation = \
                                                      stepInformationList[6]
                    currentIndex = 7
                else:
                    stepToModify.time = stepInformationList[4]
                    currentIndex = 5

                if(len(stepInformationList) >= currentIndex+1):
                    amplitudeAnswer = stepInformationList[currentIndex]
                    if(amplitudeAnswer=="flip"):
                        ## stepInformationList = [actionName, axisInfo, objectInfo, 
                        ##                        objectInformationList, timeInfo, 
                        ##                        flipAmplitudeInfo]
                        stepToModify.amplitude = "flip"
                    elif(amplitudeAnswer=="equation"):
                        ## stepInformationList = [actionName, axisInfo, objectInfo, 
                        ##                        objectInformationList, timeInfo, 
                        ##                        amplitudeTypeInfo, 
                        ##                        amplitudeEquationNameInfo]
                        stepToModify.amplitude = EquationRef()
                        stepToModify.amplitude.type = stepInformationList[currentIndex]
                        stepToModify.amplitude.equation = stepInformationList[currentIndex + 1]
                        if len(stepInformationList) >= currentIndex+3:
                            addEquation(sequence_data = sequence_data, 
                                        equationName = stepToModify.amplitude.equation)
                            sequence_data.equations[\
                                stepToModify.amplitude.equation].equation = \
                                                          stepInformationList[currentIndex + 2]
                else:
                    # print("No amplitude option added.")
                    pass

            case "rf":
                ## stepInformationList = [actionName, objectInfo, 
                ##                        objectInformationList, timeInfo, 
                ##                        addedPhaseTypeInfo, addedPhaseFloatInfo]
                stepToModify.object= stepInformationList[1]
                if stepToModify.object not in sequence_data.objects:
                    addObject(sequence_data=sequence_data, 
                              objectName=stepToModify.object,
                              typeName = "rf")
                completeObjectInformation(sequence_data = sequence_data, 
                                          objectName = stepToModify.object,
                                          objectInformationList = \
                                                         stepInformationList[2])
                
                if(len(stepInformationList) >= 7):
                    stepToModify.time = EquationRef()
                    stepToModify.time.type = stepInformationList[3]
                    stepToModify.time.equation = stepInformationList[4]
                    addEquation(sequence_data = sequence_data, 
                            equationName = stepToModify.time.equation)
                    sequence_data.equations[\
                            stepToModify.time.equation].equation = \
                                                      stepInformationList[5]
                    stepToModify.added_phase = AddedPhase()
                    stepToModify.added_phase.type = stepInformationList[6]
                    stepToModify.added_phase.float = stepInformationList[7]
                else:
                    stepToModify.time = stepInformationList[3]
                    stepToModify.added_phase = AddedPhase()
                    stepToModify.added_phase.type = stepInformationList[4]
                    stepToModify.added_phase.float = stepInformationList[5]

            case "adc":
                ## stepInformationList = [actionName, objectInfo, 
                ##                        objectInformationList, timeInfo, 
                ##                        frequencyInfo, phaseInfo,
                ##                        addedPhaseTypeInfo, addedPhaseFloatInfo,
                ##                        mdhInfoList]
                stepToModify.object= stepInformationList[1]
                if stepToModify.object not in sequence_data.objects:
                    addObject(sequence_data=sequence_data, 
                              objectName=stepToModify.object,
                              typeName = "adc")
                completeObjectInformation(sequence_data = sequence_data, 
                                          objectName = stepToModify.object,
                                          objectInformationList = \
                                                         stepInformationList[2])

                # print("+-+- len(stepInformationList) ", len(stepInformationList))
                # print("+-+- stepInformationList ", stepInformationList)
                if stepInformationList[3] == "equation":
                    stepToModify.time = EquationRef()
                    stepToModify.time.type = stepInformationList[3]
                    stepToModify.time.equation = stepInformationList[4]
                    addEquation(sequence_data = sequence_data, 
                            equationName = stepToModify.time.equation)
                    sequence_data.equations[\
                            stepToModify.time.equation].equation = \
                                                      stepInformationList[5]
                    currentIndex = 6
                else:
                    stepToModify.time = stepInformationList[3]
                    currentIndex = 4

                stepToModify.frequency = stepInformationList[currentIndex]
                stepToModify.phase= stepInformationList[currentIndex + 1]
                stepToModify.added_phase = AddedPhase()
                stepToModify.added_phase.type = stepInformationList[currentIndex + 2]
                stepToModify.added_phase.float = stepInformationList[currentIndex + 3]
                stepToModify.mdh = stepInformationList[currentIndex + 4]
            case "mark":
                ## stepInformationList = [actionName, timeInfo]
                if stepInformationList[1] == "equation":
                    stepToModify.time = EquationRef()
                    stepToModify.time.type = stepInformationList[1]
                    stepToModify.time.equation = stepInformationList[2]
                    addEquation(sequence_data = sequence_data, 
                            equationName = stepToModify.time.equation)
                    sequence_data.equations[\
                            stepToModify.time.equation].equation = \
                                                      stepInformationList[3]
                else:
                    stepToModify.time = stepInformationList[1]

            case "submit":
                pass
            case _:
                # print("The type " + actionInfo + " could not be identified.")
                pass


def completeObjectInformation(sequence_data, objectName, objectInformationList):
    """
    Completes the information for a given object and adds it to the sequence data.

    Args:
        sequence_data (SequenceData): The sequence data object.
        objectName (str): The name of the object.
        objectInformationList (list): A list containing the information for the object.

    Returns:
        None
    """
    ## objectInformationList = [typeInfo, durationInfo, objectSpecificInfo...]
    if(objectInformationList != []):
        typeInfo = objectInformationList[0]
        durationInfo = objectInformationList[1]
        match typeInfo:
            case "rf":
                ## objectInformationList = [typeInfo, durationInfo, arrayInfo, 
                ##                          arrayInformationList, initPhaseInfo, freqOffsetInfo
                ##                          thicknessInfo, flipAngleInfo, 
                ##                          purposeInfo]
                arrayInfo = objectInformationList[2]
                if arrayInfo not in sequence_data.arrays:
                    addArray(sequence_data = sequence_data, 
                             arrayName = arrayInfo)
                completeArrayInformation(sequence_data = sequence_data, 
                                         arrayName = arrayInfo,
                                         arrayInformationList = \
                                                       objectInformationList[3])
                initPhaseInfo = objectInformationList[4]
                freqOffsetInfo = objectInformationList[5]
                thicknessInfo = objectInformationList[6]
                flipangleInfo = objectInformationList[7]
                purposeInfo = objectInformationList[8]
                sequence_data.objects[objectName]=RfExcitation( 
                                                duration = durationInfo, 
                                                array = arrayInfo, 
                                                initial_phase =  initPhaseInfo,
                                                freq_offset =  freqOffsetInfo, 
                                                thickness = thicknessInfo, 
                                                flipangle = flipangleInfo, 
                                                purpose = purposeInfo) 
            case "grad":
                ## objectInformationList = [arrayInfo, durationInfo, 
                ##                          arrayName, arrayInformationList, 
                ##                          tailInfo, amplitudeInfo]
                arrayInfo = objectInformationList[2]
                if arrayInfo not in sequence_data.arrays:
                    addArray(sequence_data = sequence_data, 
                             arrayName = arrayInfo)
                completeArrayInformation(sequence_data = sequence_data, 
                                         arrayName = arrayInfo,
                                         arrayInformationList = \
                                                       objectInformationList[3])
                tailInfo = objectInformationList[4]
                amplitudeInfo = objectInformationList[5]
                sequence_data.objects[objectName]=GradientObject(
                                                    duration = durationInfo, 
                                                    array = arrayInfo,
                                                    tail = tailInfo, 
                                                    amplitude = amplitudeInfo) 
            case "adc":
                ## objectInformationList = [arrayInfo, durationInfo, samplesInfo, 
                ##                          dwelltimeInfo]
                samplesInfo = objectInformationList[2]
                dwelltimeInfo = objectInformationList[3]
                sequence_data.objects[objectName]=AdcReadout( 
                                                    duration = durationInfo, 
                                                    samples = samplesInfo, 
                                                    dwelltime = int(dwelltimeInfo)*1000) 
            case "sync":
                ## objectInformationList = [arrayInfo, durationInfo, eventInfo]
                eventInfo = objectInformationList[2]
                sequence_data.objects[objectName]=Ttl(duration = durationInfo, 
                                                      event = eventInfo) 
            case _:
                # print("The type " + typeInfo + " could not be identified.")
                pass


def completeArrayInformation(sequence_data, arrayName, arrayInformationList):
    """
    Completes the array information and adds it to the sequence_data object.

    Args:
        sequence_data (SequenceData): The sequence_data object to which the array information will be added.
        arrayName (str): The name of the array.
        arrayInformationList (list): A list containing the array information in the following format:
            [encodingInfo, typeInfo, sizeInfo, dataInfoList]

    Returns:
        None
    """
    ## arrayInformationList = [encodingInfo, typeInfo, sizeInfo, dataInfoList]
    if(arrayInformationList != []):
        encodingInfo = arrayInformationList[0]
        typeInfo = arrayInformationList[1]
        sizeInfo = int(arrayInformationList[2])
        dataInfo = arrayInformationList[3]
        sequence_data.arrays[arrayName] = GradientTemplate( encoding = encodingInfo, 
                                                            type = typeInfo, 
                                                            size = sizeInfo, 
                                                            data = dataInfo)