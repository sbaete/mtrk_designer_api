################################################################################
### mtrk project - SDL file generator designing a spin echo 2D sequence      ###
### Version 0.1.1                                                            ###
### Anais Artiges and the mtrk project team at NYU - 06/24/2025              ###
################################################################################  

from SDL_read_write.pydanticSDLHandler import *
from simpleWaveformGenerator import *
from ReadoutBlocks.mtrkReadoutBlockGenerator import add_cartesian_readout
import json
import jsbeautifier
import re

def se2d_generator():
    """
    Generates a spin echo 2D sequence data object.

    Args:

    Returns:
        SequenceData: The created sequence data object.
    """

    ### Creating the sequence data object
    sequence_data = PulseSequence()
    
    ### Filling the header part of the SDL file

    ## Sequence file and system technical specifications
    sequence_data.file = File(
        format="Custom!",
        version=1,
        measurement="abc",
        system="Skyra-XQ")

    ## Sequence-related general settings

    TE = 20000  # Echo time in us
    TR = 5000000  # Repetition time in us

    sequence_data.settings = Settings(
        readout_os=2,
        TE=TE,
        TR=TR)

    ## Information for reconstruction
    sequence_data.infos = Info(
        description="Spin Echo 2D Sequence",
        slices=1,
        fov=260,
        dz=5,
        is3D=False,
        pelines=128,
        seqstring="YARRA",
        reconstruction="%SiemensIceProgs%\\IceProgram2D")
    
    ### Calculating waveforms

    ## RF excitation and refocusing pulse

    # Generating the RF pulse
    magnitude, phase = pulse_designer("sinc", [128, 1])
    # preparing the interleaved format for SDL
    interleaved_magnitude_and_phase = [] 
    for i in range(len(magnitude)):
        interleaved_magnitude_and_phase.append(magnitude[i])
        interleaved_magnitude_and_phase.append(phase[i])

    # Generating slice selection gradient
    gamma = 42.576e6  # Hz/T
    time_bw_product = 2.7
    BW = time_bw_product/(2560e-6)  # RF bandwidth in Hz
    dz = sequence_data.infos.dz * 1e-3 # slice thickness in meters
    G_slice = BW / (gamma * dz)  # in T/m
    grad_slisel_amplitude = round(G_slice * 1e3, 2)  # convert to mT/m

    dt = 10  # 10 µs
    ramp = 400
    plateau = 256 * dt
    grad_slisel_waveform, amp, grad_slisel_ru_us, rd_us, pt_us = trap_grad(ramp, ramp, plateau, dt)

    # Generating slice refocusing gradient
    # slice_selection_half_area = (((ramp + plateau))/2)*grad_slisel_amplitude
    slice_selection_half_area = np.sum(grad_slisel_waveform)*grad_slisel_amplitude/2*10  # take discretization into account
    grad_sliref_waveform, grad_sliref_amplitude, ru_us, rd_us, pt_us = ramp_sampled_trap_grad(slice_selection_half_area, 22, 45, dt)

    # Generating crusher gradients
    dephasing = 10 * np.pi
    crusher_gradient_area = (dephasing*1e3) / (2*np.pi*gamma*1e-6 * 5e-3)  # in mT/m/s
    grad_crush_waveform, grad_crush_amplitude, ru_us, rd_us, pt_us = ramp_sampled_trap_grad(crusher_gradient_area, 22, 45, dt)
    for i, value in enumerate(grad_crush_waveform[0]):
        grad_crush_waveform[0][i] = np.round(value, 4)  # rounding to 4 decimal places for mT/m

    ### Generating sequence instructions
    
    # Synchronizing with system
    # Object
    ttl_sync = Ttl(
        duration=10,
        event="osc0")
    
    # Event
    sync_event = Sync(
        object="ttl",
        time=0)

    # Defining the gradient for slice selection excitation
    # Array
    grad_slisel = GradientTemplate(
        encoding="text",
        type="float",
        size=len(grad_slisel_waveform[0]),
        data=grad_slisel_waveform[0])

    # Object
    grad_slice_select = GradientObject(
        array="grad_slisel",
        duration=len(grad_slisel_waveform[0]) * 10,  
        tail=0,
        amplitude=grad_slisel_amplitude)
    
    # Event
    slice_gradient_event = Grad(
        axis="slice",
        object="grad_slice_select",
        time=0)
    
    # Defining the RF excitation pulse
    # Array
    rfpulse_array = GradientTemplate(
        encoding="text",
        type="complex_float",
        size=len(interleaved_magnitude_and_phase)/2,
        data=interleaved_magnitude_and_phase)  
    
    # Object
    rfpulse = RfExcitation(
        array="rfpulse",
        duration=len(interleaved_magnitude_and_phase)/2 * 20,
        initial_phase=0,
        thickness=5,
        flipangle=90,
        purpose="excitation")
    
    # Event
    rf_event = Rf(
        object="rf_excitation",
        time=(grad_slice_select.duration/2) - (rfpulse.duration/2),
        added_phase= AddedPhase(
            type="float",
            float=0))
    
    # Defining the gradient for slice selection refocusing
    # Array
    grad_sliref = GradientTemplate(
        encoding="text",
        type="float",
        size=len(grad_sliref_waveform[0]),
        data=grad_sliref_waveform[0])
    
    # Object
    grad_slice_select_refocusing = GradientObject(
        array="grad_sliref",  
        duration=len(grad_sliref_waveform[0]) * 10, 
        tail=0,
        amplitude=-grad_sliref_amplitude)

    # Event
    slice_refocus_gradient_event = Grad(
        axis="slice",
        object="grad_slice_select_refocusing",
        time=grad_slice_select.duration + 10)
    
    # Defining the gradient for slice selection refocus of se
    # Object
    grad_slice_select_refocus = GradientObject(
        array="grad_slisel",  
        duration=len(grad_slisel_waveform[0]) * 10, 
        tail=0,
        amplitude=grad_slisel_amplitude)
    
    # Event
    equation_ref_slice_ref_gradient_event = EquationRef(type = "equation",
                        equation = "equation_slice_ref_gradient_event")
    equation_slice_ref_gradient_event = Equation(equation = "set(TE)/2" + \
                        " - " + str(grad_slice_select_refocus.duration/2) + \
                        " + " + str(rfpulse.duration/2) + \
                        " + " + str(rf_event.time) + " + 10 ")
        
    slice_ref_gradient_event = Grad(
        axis="slice",
        object="grad_slice_select_refocus",
        time=equation_ref_slice_ref_gradient_event)
    
    # Defining the crusher gradients
    # Array
    grad_crush = GradientTemplate(
        encoding="text",
        type="float",
        size=len(grad_crush_waveform[0]),
        data=grad_crush_waveform[0])
    
    # Object
    grad_crusher = GradientObject(
        array="grad_crush",
        duration=len(grad_crush_waveform[0]) * 10, 
        tail=0,
        amplitude=grad_crush_amplitude)
    
    # Events
    equation_ref_crusher_gradient_event_s1 = EquationRef(type = "equation",
                        equation = "equation_crusher_gradient_event_s1")
    equation_crusher_gradient_event_s1 = Equation(equation = "set(TE)/2" + \
                        " - " + str(grad_slice_select_refocus.duration/2 + grad_crusher.duration) + \
                        " + " + str(rfpulse.duration/2) + \
                        " + " + str(rf_event.time) + " + 10 ")
    
    #start_time_crusher_gradient_event_s1 = TE/2 - (grad_slice_select_refocus.duration/2 + grad_crusher.duration) - 20
    crusher_gradient_event_s1 = Grad(
        axis="slice",
        object="grad_crusher",
        time=equation_ref_crusher_gradient_event_s1)
    crusher_gradient_event_r1 = Grad(   
        axis="read",
        object="grad_crusher",
        time=equation_ref_crusher_gradient_event_s1)
    crusher_gradient_event_p1 = Grad(   
        axis="phase",
        object="grad_crusher",
        time=equation_ref_crusher_gradient_event_s1)

    # Defining the RF refocusing pulse for se
    # Array
    rfpulse_refocusing_array = GradientTemplate(
        encoding="text",
        type="complex_float",
        size=len(interleaved_magnitude_and_phase)/2,
        data=interleaved_magnitude_and_phase)
    
    # Object
    rfpulse_refocusing = RfExcitation(
        array="rfpulse_refocusing",
        duration=len(interleaved_magnitude_and_phase)/2 * 20, 
        initial_phase=0,
        thickness=5,
        flipangle=180,
        purpose="refocusing")

    # Event 
    equation_ref_rf_ref_event= EquationRef(type = "equation",
                        equation = "equation_rf_ref_event") 
    equation_rf_ref_event = Equation(equation = "set(TE)/2" + \
                        " - " + str(rfpulse_refocusing.duration/2) + \
                        " + " + str(rfpulse.duration/2) + \
                        " + " + str(rf_event.time) + " + 10 ")
    rf_ref_event = Rf(
        object="rf_refocusing",
        time=equation_ref_rf_ref_event,
        added_phase= AddedPhase(
            type="float",
            float=0))
    
    # Defining the crusher gradients
    # Events
    equation_ref_crusher_gradient_event_s2 = EquationRef(type = "equation",
                       equation = "equation_crusher_gradient_event_s2")
    equation_crusher_gradient_event_s2 = Equation(equation = "set(TE)/2" + \
                       " + " + str(grad_slice_select_refocus.duration/2) + \
                       " + " + str(rfpulse.duration/2) + \
                       " + " + str(rf_event.time) + " + 10 ")
    crusher_gradient_event_s2 = Grad(
        axis="slice",
        object="grad_crusher",
        time=equation_ref_crusher_gradient_event_s2)
    crusher_gradient_event_r2 = Grad(
        axis="read",
        object="grad_crusher",
        time=equation_ref_crusher_gradient_event_s2)
    crusher_gradient_event_p2 = Grad(
        axis="phase",
        object="grad_crusher",
        time=equation_ref_crusher_gradient_event_s2)

    # Defining the echo time (TE) event
    equation_ref_te_event = EquationRef(type = "equation",
                        equation = "equation_te_event")
    equation_te_event = Equation(equation = "set(TE) " \
                        " + " + str(rfpulse.duration/2) + \
                        " + " + str(rf_event.time) + " + 10 ")
    te_event = Mark(time=equation_ref_te_event)
    
    
    # Defining the spoiler gradient
    # Object
    grad_spoiler = GradientObject(
        array="grad_crush",
        duration=len(grad_crush_waveform[0]) * 10, 
        tail=0,
        amplitude=grad_crush_amplitude)
    
    # Event
    # Note: the spoiler will be in a different block, the timing starts at 0.
    slice_spoiler_gradient_event = Grad(    
        axis="slice",
        object="grad_spoiler",
        time=0)

    
    # Defining the repetition time (TR) event
    #print("test " + str(TR - start_time_crusher_gradient_event_s1 - grad_crusher.duration))
    equation_ref_tr_event = EquationRef(type = "equation",
                        equation = "equation_tr_event")
    equation_tr_event = Equation(equation = "set(TR) - set(TE) " \
                        " - " + str(rfpulse.duration/2) + \
                        " - " + str(rf_event.time) + " - 10 ")
    tr_event = Mark(time=equation_ref_tr_event ) # - start_time_crusher_gradient_event_s1 - grad_crusher.duration)
    

    ## Adding objects to the sequence data
    sequence_data.objects.update({
        "rf_excitation": rfpulse,
        "rf_refocusing": rfpulse_refocusing,
        "grad_slice_select": grad_slice_select,
        "grad_slice_select_refocus": grad_slice_select_refocus,
        "grad_slice_select_refocusing": grad_slice_select_refocusing,
        "grad_crusher": grad_crusher,
        "grad_spoiler": grad_spoiler,
        "ttl": ttl_sync
    })

    ## Adding arrays to the sequence data
    sequence_data.arrays.update({
        "rfpulse":rfpulse_array,
        "rfpulse_refocusing":rfpulse_refocusing_array,
        "grad_slisel": grad_slisel,
        "grad_sliref": grad_sliref,
        "grad_crush": grad_crush
    })
    
    ## Building the sequence structure

    # Adding the main block (looping over slices) and the core blocks: 
    # block_phaseEncoding (looping over lines), and block_TR (defining the core 
    # part of the sequence)
    sequence_data.instructions.update({
        "main":{}, 
        "block_phaseEncoding":{}, 
        "block_TR":{},
        "block_SE":{},
        "block_spoiler":{}
    })
    sequence_data.instructions["main"] = Instruction(
        print_counter="on",
        print_message="Running main loop",
        steps=[])
    sequence_data.instructions["block_phaseEncoding"] = Instruction(
        print_counter="on",
        print_message="Looping over lines",
        steps=[])
    sequence_data.instructions["block_TR"] = Instruction(
        print_counter="on",
        print_message="Running TR Loop",
        steps=[])
    sequence_data.instructions["block_SE"] = Instruction(
        print_counter="on",
        print_message="Performing SE excitation and refocusing",
        steps=[])
    sequence_data.instructions["block_spoiler"] = Instruction(
        print_counter="on",
        print_message="Applying spoilers and delay to TR",
        steps=[])
    
    # Creating loop structures
    main_loop = Loop(
        counter=1,
        range=1,
        steps=[RunBlock(
            action="run_block",
            block="block_phaseEncoding")])
    sequence_data.instructions["main"].steps.append(main_loop)
    phase_encoding_loop = Loop(
        counter=2,
        range=128,
        steps=[RunBlock(block="block_TR")])
    pe_events = [
        Init(gradients="logical"),
        phase_encoding_loop,
        Submit()]
    sequence_data.instructions["block_phaseEncoding"].steps.extend(pe_events)
    spinech_events = [
        Init(gradients="logical"),
        RunBlock(block="block_SE"),
        RunBlock(block="block_spoiler"),
        Submit()]
    sequence_data.instructions["block_TR"].steps.extend(spinech_events)
   

    # Adding events to the TR block
    # List of events in the order they should be added to block_TR, this list
    # madatorily starts and ends with Init and Submit events
    se_events = [
        Init(gradients="logical"),
        sync_event,
        rf_event,
        rf_ref_event,
        slice_gradient_event,
        slice_refocus_gradient_event,
        slice_ref_gradient_event,
        crusher_gradient_event_s1,
        crusher_gradient_event_s2,
        crusher_gradient_event_r1,
        crusher_gradient_event_r2,
        crusher_gradient_event_p1,
        crusher_gradient_event_p2,
        te_event,
        Submit()
    ]
    sequence_data.instructions["block_SE"].steps.extend(se_events)

    spoiler_events = [
        Init(gradients="logical"),
        # phase_spoiler_gradient_event,
        slice_spoiler_gradient_event,
        # readout_spoiler_gradient_event,
        tr_event,
        Submit()
    ]
    sequence_data.instructions["block_spoiler"].steps.extend(spoiler_events)

    # Fillinf the equations section of the sequence data
    sequence_data.equations.update({
        "equation_crusher_gradient_event_s1": equation_crusher_gradient_event_s1,
        "equation_crusher_gradient_event_s2": equation_crusher_gradient_event_s2,
        "equation_te_event": equation_te_event,
        "equation_tr_event": equation_tr_event,
        "equation_rf_ref_event": equation_rf_ref_event,
        "equation_slice_ref_gradient_event": equation_slice_ref_gradient_event
    })

    ## Adding the cartesian readout to the prepared sequence
    add_cartesian_readout(sequence_data, 
                          "block_TR", 
                          "block_SE", 
                          "block_spoiler", 
                          sequence_data.infos.fov * 1e-2,
                          sequence_data.infos.pelines)

    return sequence_data

se2d = se2d_generator()

### writing of json schema to SDL file with formatting options
## WARNING - The path needs to be adapted to your local implementation. 
with open('se2d.mtrk', 'w') as sdlFileOut:
    options = jsbeautifier.default_options()
    options.indent_size = 4
    data_to_print = jsbeautifier.beautify(json.dumps(se2d.model_dump(mode="json")), options)
    sdlFileOut.write(re.sub(r'}, {', '},\n            {', data_to_print)) #purely aesthetic 