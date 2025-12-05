updateDynamicMeshDict() {
    # Check for required argument
    if [ -z "$1" ]; then
        echo "🛑 Error: You must provide a case directory."
        echo "Usage: updateMeshDicts <case_directory>"
        return 1
    fi

    local target_dir="$1"
    
    if [ ! -d "$target_dir" ]; then
        echo "Error: Directory '$target_dir' does not exist."
        return 1
    fi

    local FILE_COUNT=0
    
    echo "--- Starting Recursive Scan for constant/dynamicMeshDict in: $target_dir ---"

    # Find files and pipe them to the loop
    find "$target_dir" -path "*/constant/dynamicMeshDict" -type f | while IFS= read -r FILE_PATH; do
        
        FILE_COUNT=$((FILE_COUNT + 1))
        
        echo ""
        echo "--- File $FILE_COUNT Found ---"
        echo "This file will be modified: ➡️ $FILE_PATH"
        
        # FIX: Redirect the 'read' input to the terminal (/dev/tty)
        read -r -p "Do you want to apply the modifications? (Type 'y' to confirm, 'n' to skip) " response < /dev/tty

        if [[ "$response" =~ ^[Yy]$ ]]; then
            echo "Applying changes to $FILE_PATH..."
            
            # Apply the 3 sed commands
            sed -i \
                -e 's/dynamicMotionSolverListFvMesh/dynamicPointDisplacementFvMesh/g' \
                -e 's/floaterMotionCoeffs/floaterMeshMotionSolverCoeffs/g' \
                -e 's/floaterMotion;/floaterMeshMotionSolver;/g' \
                "$FILE_PATH"
            
            echo "Modification applied."
        else
            echo "Skipping $FILE_PATH."
        fi
    done

    echo ""
    echo "--- Scan Complete ---"
    echo "Finished checking $target_dir. Total files processed: $FILE_COUNT"
}