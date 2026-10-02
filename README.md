![mtrk](https://github.com/mtrk-dev/mtrk_designer_gui/blob/main/app/static/mtrk_icon.ico)

![GitHub](https://img.shields.io/github/license/mtrk-dev/mtrk_designer_api?color=%23c3ab8b)

# Description
mtrk_designer_api contains the tools to design an MRI pulse sequence, to write it as an SDL file (Json-based, .mtrk file), to convert this file to Pulseq format, and to connect to mtrk_designer_gui which is a graphical user interface (GUI) for MRI pulse sequence development. 

# Contents
This repository contains python scripts:
- miniFlashModifier (old test tool): a prototype file that reads the existing miniflash.mtrk file (available in init_data) and modifies it,
- mtrkConsoleUI: a console interface to create an SDL file,
- backendToUI: tools to connect with the GUI,
- mtrkToPulseqConverter: a tool to convert SDL files to Pulseq files using PyPulseq,
- RfPulseGenerator (WIP): a prototype to generate RF pulses,
- sdlFileCreator: SDL file generator allowing to simply define MRI pulse sequences that can be read by the mtrk project simulator and driver sequence,
- simpleWaveformGenerator: a library to generate waveforms for gradient and RF pulses,
- SDLGenerator: a file to test and call all the tools described before (main),
- se2d_generator: an example sequence with waveform calculation and SDL file creation. 

And folders:
- PrototypeFunctions: work in progress scripts to improve/extend the functionnalities of mtrk,
- SDL_read_write/pydanticSDLHandler: a set of tools to read and write SDL files using Pydantic,
- ReadoutBlocks: tools to generate readout blocks and incorporate them in existing sequence structures,
- init_data: initialization file,
- testData: example data used in the tutorial.

Additionnaly, requirements.txt helps setting the local environment by intalling the right dependencies, and Doxyfile allows to generate the doxygen documentation. 

# How to run it
An example sequence is provided in se2d_generator, py. For more details, a tutorial is available on the main project [page](https://github.com/mtrk-dev).

# Questions?
You can contact me at anais.artiges@gmail.org
