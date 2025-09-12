from mtrkReadoutBlockGenerator import *
import json
import jsbeautifier
import re

inputFilename = 'C:/Users/artiga02/mtrk_designer_gui/app/mtrk_designer_api/mtrk_designer_api.mtrk'
readoutList = ["cartesian", "radial", "spiral", "epi"]
readoutType = readoutList[0] # type of readout to add
# outputFilename = 'C:/Users/artiga02/Downloads/se2d_' + str(readoutType) + '.mtrk'
outputFilename = 'updated_sdl.mtrk'
insertion_block = "block_spinEcho" # block name to insert 
previous_block = "block_refocusing" # previous step name

## Generating a readout block and inserting it in the base sequence
## Getting fov and resolution from base sequence
def automaticReadoutBlockGenerator(readoutType = readoutType, inputFilename = inputFilename, 
                          insertion_block = insertion_block, previous_block = previous_block):
    ## Opening the base sequence file and loading it into a PulseSequence object
    with open(inputFilename) as sdlFile:
        sdlData = json.load(sdlFile)
        base_sequence = PulseSequence(**sdlData)

    ## Setting the fov and resolution according to the base sequence object
    fov = base_sequence.infos.fov * 1e-2 # imaging field of view
    resolution = base_sequence.infos.pelines # resolution
    

    ## Adding the readout block to the base sequence
    if readoutType == "cartesian":
        print("+-+-+ Generating cartesian readout")
        output_sequence = add_cartesian_readout(base_sequence, insertion_block, previous_block, fov, resolution)
    elif readoutType == "radial":
        print("+-+-+ Generating radial readout")
        output_sequence = add_radial_readout(base_sequence, insertion_block, previous_block, fov, resolution)
    elif readoutType == "spiral":
        print("+-+-+ Generating spiral readout")
        output_sequence = add_spiral_readout(base_sequence, insertion_block, previous_block, fov, resolution)
    elif readoutType == "epi":
        print("+-+-+ Generating epi readout")
        output_sequence = add_epi_readout(base_sequence, insertion_block, previous_block, fov, resolution)

    ## Generating the output sequence file
    with open(outputFilename, 'w') as sdlFileOut:
        options = jsbeautifier.default_options()
        options.indent_size = 4
        data_to_print = jsbeautifier.beautify(json.dumps(output_sequence.model_dump(mode="json")), options)
        sdlFileOut.write(re.sub(r'}, {', '},\n            {', data_to_print)) #purely aesthetic 

## Generating a readout block and inserting it in the base sequence
## Force setting fov and resolution
def manualReadoutBlockGenerator(readoutType = readoutType, inputFilename = inputFilename, 
                          insertion_block = insertion_block, previous_block = previous_block,
                          fov = 260, resolution = 128):
    ## Opening the base sequence file and loading it into a PulseSequence object
    with open(inputFilename) as sdlFile:
        sdlData = json.load(sdlFile)
        base_sequence = PulseSequence(**sdlData)

    ## Setting the base sequence fov and resolution according to the provided values
    base_sequence.infos.fov = fov # imaging field of view
    fov = fov * 1e-2
    base_sequence.infos.pelines = resolution # resolution (warning, this is not changing the iteration number for now)
    
    ## Adding the readout block to the base sequence
    if readoutType == "cartesian":
        output_sequence = add_cartesian_readout(base_sequence, insertion_block, previous_block, fov, resolution)
    elif readoutType == "radial":
        output_sequence = add_radial_readout(base_sequence, insertion_block, previous_block, fov, resolution)
    elif readoutType == "spiral":
        output_sequence = add_spiral_readout(base_sequence, insertion_block, previous_block, fov, resolution)
    elif readoutType == "epi":
        output_sequence = add_epi_readout(base_sequence, insertion_block, previous_block, fov, resolution)

    ## Generating the output sequence file
    with open(outputFilename, 'w') as sdlFileOut:
        options = jsbeautifier.default_options()
        options.indent_size = 4
        data_to_print = jsbeautifier.beautify(json.dumps(output_sequence.model_dump(mode="json")), options)
        sdlFileOut.write(re.sub(r'}, {', '},\n            {', data_to_print)) #purely aesthetic

## Testing functions

## Defining parameters for the test
## Res 128 TE200 TR4000 single line
# inputFilename = 'C:/Users/artiga02/mtrk_designer_gui/app/mtrk_designer_api/base_sequence_singleLine_TE200_TR4000_res128.mtrk'
## Res 128 TE60 TR4000 single line
# inputFilename = 'C:/Users/artiga02/mtrk_designer_gui/app/mtrk_designer_api/base_sequence_singleLine_TE60_TR4000_res128.mtrk'
## Res 64 TE60 TR4000 single line
# inputFilename = 'C:/Users/artiga02/mtrk_designer_gui/app/mtrk_designer_api/base_sequence_singleLine_TE60_TR4000_res64.mtrk'
## Res 128 TE200 TR4000 multi lines
# inputFilename = 'C:/Users/artiga02/mtrk_designer_gui/app/mtrk_designer_api/base_sequence_multiLines_TE200_TR4000_res128.mtrk'
## Res 128 TE60 TR4000 multi lines
#inputFilename = 'C:/Users/artiga02/ISMRM-2025-Surfing-School-Hands-On-Open-Source-MR/Sequences/base_sequences/base_sequence_singleShot_TE60_TR4000_res128.mtrk'
## Res 64 TE60 TR4000 multi lines
# inputFilename = 'C:/Users/artiga02/mtrk_designer_gui/app/mtrk_designer_api/base_sequence_multiLines_TE60_TR4000_res64.mtrk'

# readoutList = ["cartesian", "radial", "spiral", "epi"]
# readoutType = readoutList[2] # type of readout to add
# outputFilename = 'C:/Users/artiga02/ISMRM-2025-Surfing-School-Hands-On-Open-Source-MR/Sequences/mtrk/se2d_' + str(readoutType) + '_shortTE.mtrk'
# insertion_block = "block_spinEcho" # block name to insert
# previous_block = "block_refocusing" # previous step name

# ## Testing the automatic readout block generator
#automaticReadoutBlockGenerator(readoutType = readoutType, inputFilename = inputFilename,
                         #insertion_block = insertion_block, previous_block = previous_block)

## Defining extra parameters for the test
# fov = 260 # imaging field of view
# resolution = 128 # resolution (warning, this is not changing the iteration number for now)

# ## Testing the manual readout block generator
# manualReadoutBlockGenerator(readoutType = readoutType, inputFilename = inputFilename,
#                          insertion_block = insertion_block, previous_block = previous_block,
#                          fov = fov, resolution = resolution)
