import fileinput, math, re
import numpy as np

x1 = 80
x2 = 120
y1 = 30
y2 = 75
z1 = 95
z2 = 105

dxA = 40
dyA = 45
dzA = -20

rAp1 = np.array([x1,y1,z1])
a1 = np.array([x1-dxA , y1-dyA, z1+dzA])
rAp2 = np.array([x2,y1,z1])
a2 = np.array([x2+dxA, y1-dyA, z1+dzA])
rAp3 = np.array([x2,y2,z1])
a3 = np.array([x2+dxA, y2+dyA, z1+dzA])
rAp4 = np.array([x1,y2,z1])
a4 = np.array([x1-dxA, y2+dyA, z1+dzA])

rl1 = np.linalg.norm(rAp1-a1)
rl2 = np.linalg.norm(rAp2-a2)
rl3 = np.linalg.norm(rAp3-a3)
rl4 = np.linalg.norm(rAp4-a4)

#Making attachment points relative to body centre
def get_centre_of_rotation(file_path):

    # Initialize an empty list to store the line containing "centreOfRotation"
    centre_of_rotation_line = None

    # Open the file and read line by line
    with open(file_path, 'r') as file:
        for line in file:
            if 'centreOfRotation' in line:
                centre_of_rotation_line = line.strip()
                break

    # Check if the line was found
    if centre_of_rotation_line:
        # Extract the numbers from the line
        start_idx = centre_of_rotation_line.find('(') + 1
        end_idx = centre_of_rotation_line.find(')')
        numbers_str = centre_of_rotation_line[start_idx:end_idx].strip()
        
        # Convert the string numbers to a list of floats
        numbers_list = [float(num) for num in numbers_str.split()]
        
        # Convert the list to a numpy array
        centre_of_rotation_array = np.array(numbers_list)
        
        return centre_of_rotation_array
    else:
        raise ValueError("The line containing 'centreOfRotation' was not found.")

CoR = get_centre_of_rotation('0/uniform/floaterMotionState')
rAp1 = rAp1 - CoR
rAp2 = rAp2 - CoR
rAp3 = rAp3 - CoR
rAp4 = rAp4 - CoR

with open('constant/mooringDict', 'w') as moorFil:
    moorFil.write('mooringLine1\n');
    moorFil.write('{\n');
    moorFil.write('    floaterMotionRestraint linearSpring;\n');
    moorFil.write('    anchor          (%g %g %g);\n' % (a1[0], a1[1], a1[2]));
    moorFil.write('    refAttachmentPt (%g %g %g);\n' % (rAp1[0], rAp1[1], rAp1[2]));
    moorFil.write('    restLength      %g;\n' % rl1)

    moorFil.write('    #include "moorParms"\n');
    moorFil.write('}\n\n');

    moorFil.write('mooringLine2\n');
    moorFil.write('{\n');
    moorFil.write('    $mooringLine1;\n');
    moorFil.write('    anchor          (%g %g %g);\n' % (a2[0], a2[1], a2[2]));
    moorFil.write('    refAttachmentPt (%g %g %g);\n' % (rAp2[0], rAp2[1], rAp2[2]));
    moorFil.write('    restLength      %g;\n' % rl2)
    moorFil.write('}\n\n');

    moorFil.write('mooringLine3\n');
    moorFil.write('{\n');
    moorFil.write('    $mooringLine1;\n');
    moorFil.write('    anchor          (%g %g %g);\n' % (a3[0], a3[1], a3[2]));
    moorFil.write('    refAttachmentPt (%g %g %g);\n' % (rAp3[0], rAp3[1], rAp3[2]));
    moorFil.write('    restLength      %g;\n' % rl3)
    moorFil.write('}\n\n');

    moorFil.write('mooringLine4\n');
    moorFil.write('{\n');
    moorFil.write('    $mooringLine1;\n');
    moorFil.write('    anchor          (%g %g %g);\n' % (a4[0], a4[1], a4[2]));
    moorFil.write('    refAttachmentPt (%g %g %g);\n' % (rAp4[0], rAp4[1], rAp4[2]));
    moorFil.write('    restLength      %g;\n' % rl4)
    moorFil.write('}\n\n');
    """
    moorFil.write('mooringLine5\n');
    moorFil.write('{\n');
    moorFil.write('    $mooringLine1;\n');
    moorFil.write('    anchor          (%g %g %g);\n' % (a5[0], a5[1], a5[2]));
    moorFil.write('    refAttachmentPt (%g %g %g);\n' % (m5[0], m5[1], m5[2]));
    moorFil.write('}\n\n');

    moorFil.write('mooringLine6\n');
    moorFil.write('{\n');
    moorFil.write('    $mooringLine1;\n');
    moorFil.write('    anchor          (%g %g %g);\n' % (a6[0], a6[1], a6[2]));
    moorFil.write('    refAttachmentPt (%g %g %g);\n' % (m6[0], m6[1], m6[2]));
    moorFil.write('}\n\n');

    moorFil.write('mooringLine7\n');
    moorFil.write('{\n');
    moorFil.write('    $mooringLine1;\n');
    moorFil.write('    anchor          (%g %g %g);\n' % (a7[0], a7[1], a7[2]));
    moorFil.write('    refAttachmentPt (%g %g %g);\n' % (m7[0], m7[1], m7[2]));
    moorFil.write('}\n\n');

    moorFil.write('mooringLine8\n');
    moorFil.write('{\n');
    moorFil.write('    $mooringLine1;\n');
    moorFil.write('    anchor          (%g %g %g);\n' % (a8[0], a8[1], a8[2]));
    moorFil.write('    refAttachmentPt (%g %g %g);\n' % (m8[0], m8[1], m8[2]));
    moorFil.write('}');
    """
