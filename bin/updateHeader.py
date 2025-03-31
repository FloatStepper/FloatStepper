#!/usr/bin/env python3
import os
import re

def replace_header(file_path, new_header):
    """
    Replaces the header in a file with the new header.
    The header is identified by the first occurrence of text between
    /*--------------------------------------------------- -----------------------*\
    and 
    \*--------------------------------------------------- -----------------------*/
    """
    with open(file_path, 'r') as f:
        content = f.read()
    
    # Pattern to match header (with any content between the markers)
    header_pattern = r'/\*---------------------------------------------------------------------------\*\\(.*?)\\*---------------------------------------------------------------------------\*/'
    # Use DOTALL flag to make . match newlines as well
    match = re.search(header_pattern, content, re.DOTALL)
    
    if match:
        # Replace the header with the new header
        modified_content = content.replace(match.group(0), new_header)
        
        # Write the modified content back to the file
        with open(file_path, 'w') as f:
            f.write(modified_content)
        
        print(f"Header replaced in {file_path}")
    else:
        # If no header is found, do nothing
        print(f"No header found in {file_path}")

def main():
    # Read the new header content
    try:
        with open("header", 'r') as f:
            new_header = f.read().strip()
    except FileNotFoundError:
        print("Error: header file not found")
        return
    
    # Walk through all subdirectories
    for root, dirs, files in os.walk('..'):
        for file in files:
            # Check if file has .C or .H extension (case sensitive)
            if file.endswith('.C') or file.endswith('.H'):
                file_path = os.path.join(root, file)
                replace_header(file_path, new_header)

if __name__ == "__main__":
    main()