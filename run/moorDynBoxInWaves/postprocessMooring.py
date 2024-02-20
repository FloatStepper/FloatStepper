import os
import numpy as np
import datetime

# Define the parameters
Filelocation = 'Mooring/Lines1.out'  # Replace with the actual file name
numLines = 4
nline = numLines
model = 'Catenary Mooring SPAR case'
simdate = datetime.date.today().strftime("%d-%b-%Y")

# Create the output directory
output_directory = "ParaviewData"
os.makedirs(output_directory, exist_ok=True)

# Read the data file
data = np.genfromtxt(Filelocation, skip_header=2, delimiter='\t')

# Extract time values and node coordinates
time_values = data[:, 0]
node_coordinates = data[:, 1:10].reshape(-1, 3)  # Assuming each node has 3 coordinates

# Create the output directory for VTK files
vtk_output_directory = os.path.join(output_directory, "VTKFiles")
os.makedirs(vtk_output_directory, exist_ok=True)

# Process the data for each time step
for it, current_time in enumerate(time_values):
    output_filename = os.path.join(vtk_output_directory, f'mooring_{it}.vtp')

    with open(output_filename, 'w') as f:
        f.write('<?xml version="1.0"?>\n')
        f.write('<!-- foamStar Visualization using ParaView -->\n')
        f.write(f'<!--   model: {model} - ran on {simdate} -->\n')
        f.write('<!--   mooring:  MoorDyn -->\n')
        f.write(f'<!--   time:  {current_time:.7f} -->\n')
        f.write('<VTKFile type="PolyData" version="0.1">\n')
        f.write('  <PolyData>\n')
        f.write(f'    <Piece NumberOfPoints="{len(node_coordinates)}" NumberOfLines="{len(node_coordinates) - 1}">\n')
        f.write('      <Points>\n')
        f.write('        <DataArray type="Float32" NumberOfComponents="3" format="ascii">\n')

        for pt in node_coordinates:
            f.write(f'          {pt[0]:.5f} {pt[1]:.5f} {pt[2]:.5f}\n')

        f.write('        </DataArray>\n')
        f.write('      </Points>\n')
        f.write('      <Lines>\n')
        f.write('        <DataArray type="Int32" Name="connectivity" format="ascii">\n')

        for i in range(len(node_coordinates) - 1):
            f.write(f'          {i} {i + 1}\n')

        f.write('        </DataArray>\n')
        f.write('        <DataArray type="Int32" Name="offsets" format="ascii">\n')

        for i in range(len(node_coordinates) - 1):
            n = 2 * i
            f.write(f'          {n}\n')

        f.write('        </DataArray>\n')
        f.write('      </Lines>\n')
        f.write('    </Piece>\n')
        f.write('  </PolyData>\n')
        f.write('</VTKFile>\n')

    print(f'Processed Time Step {it}, Output File: {output_filename}')

# Create the PVD file
pvd_filename = os.path.join(output_directory, 'Mooring.pvd')
with open(pvd_filename, 'w') as pvd_file:
    pvd_file.write('<?xml version="1.0"?>\n')
    pvd_file.write('<VTKFile type="Collection" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">\n')
    pvd_file.write('  <Collection>\n')

    for it in range(len(time_values)):
        pvd_file.write(f'    <DataSet timestep="{time_values[it]:.7f}" group="" part="0"\n')
        pvd_file.write(f'             file="VTKFiles/mooring_{it}.vtp"/>\n')

    pvd_file.write('  </Collection>\n')
    pvd_file.write('</VTKFile>\n')

print("Conversion completed successfully.")
