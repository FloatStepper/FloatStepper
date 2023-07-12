import os

def add_header_to_files(directory, header_file):
    # Get the content of the header file
    with open(header_file, 'r') as f:
        header_content = f.read()

    # Iterate over all files in the directory and its subdirectories
    for root, dirs, files in os.walk(directory):
        for file in files:
            # Check if the file is a source file (C or H)
            if file.endswith('.C') or file.endswith('.H'):
                file_path = os.path.join(root, file)

                # Read the existing content of the file
                with open(file_path, 'r') as f:
                    file_content = f.read()

                # Check if the header is already present in the file
                if header_content in file_content:
                    print(f'Header already exists in file: {file_path}')
                else:
                    # Combine the header and file content
                    new_content = header_content + '\n' + file_content

                    # Write the new content back to the file
                    with open(file_path, 'w') as f:
                        f.write(new_content)

                    print(f'Header added to file: {file_path}')

# Usage example
directory = '../src'  # Replace with your directory path
header_file = 'header'  # Replace with your header file path
add_header_to_files(directory, header_file)
directory = '../applications'  # Replace with your directory path
add_header_to_files(directory, header_file)
