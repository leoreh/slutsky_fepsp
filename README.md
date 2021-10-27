# Welcome to the Slutsky Lab fEPSP Repository!

Here you will find resources for recording *in vivo* field excitatory post synaptic potentials (fEPSPs) and a MATLAB© based package for analysis of the recorded signals.

## Resources

- Stimulation protocols written for winWCP and Arduino
- A spreadsheet for easy localization of implant coordinates during surgery
- CAD files for printing and/or machining components designed for custom head fixation apparatus



## Analysis

**See the ExampleData directory for a full rundown of the analysis functions.**

#### Pipeline

**Before you begin** - Prepare the raw data including times of stimulus onset.<br>

<img src="/Graphics/rawData.png" width=100 ><br>

![](Graphics/rawData.png)

**fepsp_org2traces** - Organize the raw data according to the stimulus protocol, number of channels, number of intensities, etc..<br>
format accepted by Take a raw channel data, extract only the stimulus protocol repetitions, and organize everything in a way the package expect.<br>

<img src="Graphics/traces.tif" width=100 ><br>

**fepsp_markings** - Plot traces and manually mark the start and peak of the evoked responses. <br>

<img src="Graphics/markings.tif" width=100 ><br>

**fepsp_analyse** - Calculate slope and amplitude for each trace and the average fepsp_org2traces.<br>
**fepsp_summaryPlot** - Summary plot of the data after analysis.

<img src="Graphics/summaryPlot.tif" width=100 ><br>


#### Requirements

- Hardware: This package can run on any computer capable of running MATLAB 2019a.
- Software: This package requires MATLAB© 2019a or higher with the following toolboxes:
 - Image Processing Toolbox™,
 - Signal Processing Toolbox™
 - Statistics and Machine Learning Toolbox™

#### Citation
The analysis package was written by Lior de Marcas and Leore R. Heim.<br>
Cite as: Heim et al., STAR protocol.
