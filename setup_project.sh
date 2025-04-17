#!/bin/bash
# File: nano_analysis/setup_project.sh

# Help function
show_help() {
    echo "Usage: ./setup_project.sh [OPTIONS]"
    echo ""
    echo "Setup a new Nanopore analysis project."
    echo ""
    echo "Options:"
    echo "  -n, --name NAME      Project name (required)"
    echo "  -r, --reference PATH Path to reference genome (required)"
    echo "  -d, --data PATH      Path to data directory with barcode folders"
    echo "  -h, --help           Show this help message"
    echo ""
    echo "Example:"
    echo "  ./setup_project.sh -n my_project -r /path/to/reference.fasta -d /path/to/nanopore/output"
}

# Parse command line arguments
PROJECT_NAME=""
REFERENCE_PATH=""
DATA_PATH=""

while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -n|--name)
            PROJECT_NAME="$2"
            shift
            shift
            ;;
        -r|--reference)
            REFERENCE_PATH="$2"
            shift
            shift
            ;;
        -d|--data)
            DATA_PATH="$2"
            shift
            shift
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            echo "Error: Unknown option $1"
            show_help
            exit 1
            ;;
    esac
done

# Check required parameters
if [ -z "$PROJECT_NAME" ] || [ -z "$REFERENCE_PATH" ]; then
    echo "Error: Project name and reference path are required."
    show_help
    exit 1
fi

if [ ! -f "$REFERENCE_PATH" ]; then
    echo "Error: Reference file does not exist: $REFERENCE_PATH"
    exit 1
fi

# Create project directory structure
PROJECT_DIR="projects/$PROJECT_NAME"
echo "Creating project directory structure in $PROJECT_DIR..."

mkdir -p "$PROJECT_DIR/data"
mkdir -p "$PROJECT_DIR/reference"
mkdir -p "$PROJECT_DIR/results"

# Copy reference genome
REFERENCE_FILENAME=$(basename "$REFERENCE_PATH")
echo "Copying reference genome to $PROJECT_DIR/reference/$REFERENCE_FILENAME..."
cp "$REFERENCE_PATH" "$PROJECT_DIR/reference/"

# Create project config
CONFIG_FILE="$PROJECT_DIR/config.yaml"
echo "Creating project configuration file: $CONFIG_FILE..."

cat > "$CONFIG_FILE" << EOL
# Project information
project_name: "$PROJECT_NAME"

# Input/Output paths
input_dir: "projects/$PROJECT_NAME/data"
reference_genome: "projects/$PROJECT_NAME/reference/$REFERENCE_FILENAME"
reference_name: "reference"
results_dir: "projects/$PROJECT_NAME/results"

# Processing parameters
threads: 16
min_read_quality: 8
min_read_length: 500
headcrop: 50
variant_depth: 30
variant_quality: 60
EOL

# Link data if provided
if [ ! -z "$DATA_PATH" ]; then
    if [ -d "$DATA_PATH" ]; then
        echo "Creating symbolic links to data from $DATA_PATH..."
        for dir in "$DATA_PATH"/barcode*; do
            if [ -d "$dir" ]; then
                ln -sf "$(realpath "$dir")" "$PROJECT_DIR/data/"
            fi
        done
    else
        echo "Warning: Data directory does not exist: $DATA_PATH"
        echo "You'll need to manually add your data to $PROJECT_DIR/data/"
    fi
fi

echo ""
echo "Project setup complete!"
echo ""
echo "To run the analysis:"
echo "  snakemake --configfile $CONFIG_FILE --use-conda --cores 16 --latency-wait 60"
echo ""
echo "To customize parameters, edit $CONFIG_FILE"