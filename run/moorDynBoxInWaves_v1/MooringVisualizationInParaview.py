# This code to run to be placed in main folder where it can see 'Mooring/' folder
import os
import numpy as np
import time
from datetime import datetime

# Input Parameters
numLines = 4 # Number of mooring lines 
ParaviewStart = 0 # Paraview Visualiasation starting time
ParaviewEnd = 10 # End time of Simulation you needed for Visualisation
dtParaview = 0.1 # Timestep of time folders in the visualisation
Scale = 1 # Unless you have scaled the mooring lines 

# Start timer for tracking execution time
current_dir = os.path.dirname(os.path.realpath(__file__))  # Get current directory
mooring_dir = os.path.join(current_dir, 'Mooring')  # Navigate into 'mooring' folder
start_time = time.time()
nline = numLines
model = 'Developed by Sithik Aliyar - Mooring Visualisation'
simdate =  datetime.now().date()
ParaviewdataStoring = mooring_dir  # Store data in 'Mooring' folder



print("Reading MoorDyn *.out files...")

# Generating number of nodes and number of segments
obj = {'moorDyn': {}}  # Create an empty dictionary to store mooring data

for iline in range(1, nline + 1):
    obj['moorDyn']['Line{}'.format(iline)] = {}
    filename = os.path.join(mooring_dir, 'Line{}.out'.format(iline))
    
    try:
        with open(filename, 'r') as fid:
            header = fid.readline().strip().split()
            data = np.loadtxt(fid, skiprows=2)
            
            # Store the data in the structure
            for icol, name in enumerate(header):
                obj['moorDyn']['Line{}'.format(iline)][name] = data[:, icol]
        print(f"Successfully read Line{iline}.out")
    except IOError:
        print(f"\nNo MoorDyn *.out file saved for Line{iline}")

Moordyn_time_vector = data[:, 0]

print(f"Finished reading MoorDyn files. Time taken: {time.time() - start_time:.2f} seconds.")

# Generating nnodes and nsegments from obj struct
nnode = []
nsegment = []

for k in obj['moorDyn'].keys():
    Fields = obj['moorDyn'][k]
    Linefieldnames = Fields.keys()
    OverallLength = len(Linefieldnames) - 1
    nnode.append((OverallLength + 1) // 4)
    nsegment.append(nnode[-1] - 1)

# Writing in VTP format
Paraviewtimes = np.arange(ParaviewStart, ParaviewEnd + dtParaview, dtParaview)

vtpindex = []
for it, time_val in enumerate(Paraviewtimes):
    index = np.argmin(np.abs(Moordyn_time_vector - time_val))
    vtpindex.append(index)

print("Starting to write VTP files...")

for it, ind in enumerate(vtpindex):
    current_time = Paraviewtimes[it]
    filename = os.path.join(ParaviewdataStoring, 'mooring_{}.vtp'.format(ind))
    
    with open(filename, 'w') as fid:
        fid.write('<?xml version="1.0"?>\n')
        fid.write('<!-- FloatStepper Visualization using ParaView -->\n')
        fid.write(f'<!--   model: {model} - ran on {simdate} -->\n')
        fid.write('<!--   mooring:  MoorDyn -->\n')
        fid.write(f'<!--   time:  {Moordyn_time_vector[ind]} -->\n')
        fid.write('<VTKFile type="PolyData" version="0.1">\n')
        fid.write('  <PolyData>\n')

        for iline in range(1, nline + 1):
            fid.write(f'    <Piece NumberOfPoints="{nnode[iline - 1]}" NumberOfLines="{nsegment[iline - 1]}">\n')

            # Write points
            fid.write('      <Points>\n')
            fid.write('        <DataArray type="Float32" NumberOfComponents="3" format="ascii">\n')
            for inode in range(nnode[iline - 1]):
                pt = [obj['moorDyn']['Line{}'.format(iline)]['Node{}px'.format(inode)][ind],
                      obj['moorDyn']['Line{}'.format(iline)]['Node{}py'.format(inode)][ind],
                      obj['moorDyn']['Line{}'.format(iline)]['Node{}pz'.format(inode)][ind]]
                fid.write(f'          {pt[0]:.5f} {pt[1]:.5f} {pt[2]:.5f}\n')
            fid.write('        </DataArray>\n')
            fid.write('      </Points>\n')

            # Write lines connectivity
            fid.write('      <Lines>\n')
            fid.write('        <DataArray type="Int32" Name="connectivity" format="ascii">\n')
            count = 0
            for isegment in range(1, nsegment[iline - 1] + 1):
                fid.write(f'          {count} {count + 1}\n')
                count += 1
            fid.write('        </DataArray>\n')

            # Write offsets
            fid.write('        <DataArray type="Int32" Name="offsets" format="ascii">\n')
            fid.write('         ')
            for isegment in range(1, nsegment[iline - 1] + 1):
                n = 2 * isegment
                fid.write(f' {n}')
            fid.write('\n')
            fid.write('        </DataArray>\n')
            fid.write('      </Lines>\n')

            # Write cell data
            fid.write('      <CellData>\n')
            fid.write('        <DataArray type="Float32" Name="Segment Tension" NumberOfComponents="1" format="ascii">\n')
            for isegment in range(nsegment[iline - 1]):
                fid.write(f'          {obj["moorDyn"]["Line{}".format(iline)]["Seg{}Te".format(isegment)][ind]}\n')
            fid.write('        </DataArray>\n')
            fid.write('      </CellData>\n')

            fid.write('    </Piece>\n')

        fid.write('  </PolyData>\n')
        fid.write('</VTKFile>\n')

    # Show progress every 5 files
    if (it + 1) % 5 == 0 or it == len(vtpindex) - 1:
        print(f"Written {it + 1} VTP files out of {len(vtpindex)}. Time: {time.time() - start_time:.2f} seconds.")

print(f"Finished writing all VTP files. Total time: {time.time() - start_time:.2f} seconds.")

# To create a pvd file to view the Mooring lines as matching with timestep
filename = os.path.join(ParaviewdataStoring, 'Mooring.pvd')

with open(filename, 'w') as fid:
    fid.write('<?xml version="1.0"?>\n')
    fid.write(' <VTKFile type="Collection" version="0.1" \n')
    fid.write(' byte_order="LittleEndian" \n')
    fid.write(' compressor="vtkZLibDataCompressor"> \n')
    fid.write('  <Collection>\n')

    for it, current_time in enumerate(Paraviewtimes):
        index = np.argmin(np.abs(data[:, 0] - current_time))
        MatchingTimestep = index
        fid.write(f'<DataSet timestep="{current_time}" group="" part="0" \n')
        filename1 = os.path.join(".", f'mooring_{MatchingTimestep}')
        fid.write(f'file="{filename1}.vtp"/>\n')

    fid.write('  </Collection>\n')
    fid.write('</VTKFile>\n')

print("Finished writing PVD file.")
