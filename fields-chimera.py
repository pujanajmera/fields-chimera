import numpy as np
import argparse
import warnings

def process_data(file_path):
    '''
    This function processes all data from the field file
    Inputs:
        file_path: Path to the field file
    Outputs:
        sample_density: Sample density
        volume_box: Volume of the box in Angstroms
        center: Center of the box in Angstroms
        basis_matrix: Basis matrix used for rotation
        field: Electric field values
    '''
    # Initialize variables
    sample_density = []
    volume_box = []
    center = []
    basis_matrix = []
    field = []

    # Read the file
    reading_basis_matrix = False
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('#Sample Density'):
                parts = line.split(';')
                try:
                    sample_density = [float(x) for x in parts[0].split()[2:5]]
                except:
                    raise ValueError("Sample density line not formatted correctly")
                try:
                    volume_box = [float(x) for x in parts[1].split()[2:5]]
                except:
                    raise ValueError("Volume box part of density line not formatted correctly")
            elif line.startswith('#Center'):
                try:
                    center = [float(x) for x in line.split()[1:4]]
                except:
                    raise ValueError("Center line not formatted correctly")
            elif line.startswith('#Basis Matrix'):
                reading_basis_matrix = True
            elif reading_basis_matrix and line.startswith('#'):
                try:
                    basis_matrix.append([float(x) for x in line[1:].split()[0:3]])
                except:
                    raise ValueError(f"Basis matrix line not formatted correctly; matrix list of length {len(basis_matrix)}")
                if len(basis_matrix) == 3:
                    reading_basis_matrix = False
            elif not line.startswith('#'):
                try:
                    field.append([float(x) for x in line.split()])
                except:
                    raise ValueError(f"Field line not formatted correctly; field list of length {len(field)}")
    if not sample_density:
        warnings.warn("Sample density not found, ignoring for checking")
    if not field:
        raise ValueError("No field data found at all, exiting...")
    if not center or not basis_matrix:
        warnings.warn("Center or basis matrix not found, no transformation will be made")
    if sample_density:
        check_field(np.array(field), np.array(sample_density))
    print(np.array(basis_matrix))
    return np.array(sample_density), np.array(volume_box), np.array(center), np.array(basis_matrix), np.array(field)

def check_field(field_array, sample_density_array):
    '''
    This function checks if the field provided matches the sample density
    Inputs:
        field_array: Field array
        sample_density_array: Sample density array
    Outputs:
        None
    '''
    print(sample_density_array)
    print(field_array.shape)
    if field_array.shape[0] != np.product(np.concatenate([(2*sample_density_array+1)[0:2],np.expand_dims(2*sample_density_array[2]+1, axis=0)])):
        raise ValueError(f"Field provided does not match sample density, field of shape {field_array.shape[0]} does not match expected sample amount of {np.product(np.concatenate([(2*sample_density_array+1)[0:2],np.expand_dims(sample_density_array[2], axis=0)]))}")
    else:
        print("Field matches sample density, check passed, continuing...")

def transform_field(field_array, center_array, basis_matrix_array):
    '''
    This function transforms the field to the new basis and new center
    Inputs:
        field_array: Field array
        center_array: Center array
        basis_matrix_array: Basis matrix array
    Outputs:
        transformed_field_array: Transformed field array
    '''
    field_coords = field_array[:, 0:3]
    field_vecs = field_array[:, 3:6]
    print(field_coords, field_vecs, basis_matrix_array)
    # Need to transform and recenter field coords, but only transform field vecs
    transformed_field_coords = field_coords@basis_matrix_array.T
    transformed_field_coords = transformed_field_coords + center_array
    transformed_field_vecs = np.matmul(field_vecs, basis_matrix_array.T)
    transformed_field_array = np.concatenate((transformed_field_coords, transformed_field_vecs), axis=1)
    return transformed_field_array

def generate_tip_tail_vectors(field_array, sample_density_array, volume_box_array):
    '''
    This function generates the tip and tail vectors for the field
    Inputs:
        field_array: Field array (transformed or not) Shape (N, 6)
        sample_density_array: Sample density array Shape (3,)
        volume_box_array: Volume box array
    Outputs:
        tip_to_tail_vectors: Tip to tail vectors
    '''
    if sample_density_array!=[] and volume_box_array!=[]:
        max_vector_length = 0.9*np.min(volume_box_array/sample_density_array)
    else:
        max_vector_length = 0.25

    # Normalize the vectors in field_array
    field_array[:, 3:6] = field_array[:, 3:6]/np.linalg.norm(field_array[:, 3:6], axis=1)[:, None] * max_vector_length
    tail_vecs = field_array[:, 0:3] - field_array[:, 3:6]/2
    tip_vecs = field_array[:, 0:3] + field_array[:, 3:6]/2

    tip_to_tail_vectors = np.concatenate((tail_vecs, tip_vecs), axis=1)

    return tip_to_tail_vectors

def generate_bild_file(tip_to_tail_vectors, 
                       transformed_field_array, 
                       chimera_type, 
                       percentile, 
                       sparsify_factor,
                       output):
    '''
    This function generates the BILD file for the 3D vector field
    Inputs:
        tip_to_tail_vectors: Tip to tail vectors
        transformed_field_array: Field array
        chimera_type: Chimera type
        percentile: Percentile cutoff for sparsifying the field
        sparsify_factor: Sparsification factor for the field
    Outputs:
        None
    '''
    transformed_field_vecs = transformed_field_array[:, 3:6]
    field_mags = np.linalg.norm(transformed_field_vecs, axis=1)
    field_mags_max = np.max(field_mags)
    field_mags_min = np.min(field_mags)
    field_mags = (field_mags-field_mags_min)/field_mags_max
    percentile_cutoff = np.percentile(field_mags, percentile)
    r = (1-field_mags)[::sparsify_factor]
    g = (1-field_mags)[::sparsify_factor]
    b = (np.ones((field_mags.shape[0], 1))).astype(int)[::sparsify_factor]
    field_mags = field_mags[::sparsify_factor]
    tip_to_tail_vectors = tip_to_tail_vectors[::sparsify_factor]
    with open(f'{output}.bild', 'w') as bild:
        bild.write(".transparency 0.25\n")
        for i in range(tip_to_tail_vectors.shape[0]):
            if field_mags[i] > percentile_cutoff:
                bild.write(f".color {r[i]} {g[i]} {b[i][0]}\n")
                bild.write(f".arrow {tip_to_tail_vectors[i, 0]} {tip_to_tail_vectors[i, 1]} {tip_to_tail_vectors[i, 2]} {tip_to_tail_vectors[i, 3]} {tip_to_tail_vectors[i, 4]} {tip_to_tail_vectors[i, 5]} 0.01 0.02 0.5\n")
    bild.close()

def main():
    parser = argparse.ArgumentParser(description='Process field data')
    parser.add_argument('-v', metavar='field_file_path', type=str, help="Input vector field file")
    parser.add_argument('-p', metavar='protein_file_path', type=str, help="Input protein PDB file")
    parser.add_argument('-o', metavar='output_file_name', type=str, help="Output file name (will end in .bild by default, so don't include an extension)", default='output')
    parser.add_argument('-ch', metavar='chimera_type', type=str, help="Output file type (either 'chimera' or 'chimerax'); Default: chimerax", default='chimerax')
    parser.add_argument('-c', metavar='percentile', type=int, help="Top Nth percentile for sparsifying the field; Default: 0", default=0)
    parser.add_argument('-s', metavar='sparsify_factor', type=int, help="Sparsification factor for the field (can be used in tandem with percentile_cutoff); Default: 1", default=1)

    args = parser.parse_args()
    sample_density_array, volume_box_array, center_array, basis_matrix_array, field_array = process_data(args.v)

    if center_array!=[] and basis_matrix_array!=[]:
        transformed_field_array = transform_field(field_array, center_array, basis_matrix_array)
    else:
        transformed_field_array = field_array
    
    tip_to_tail_vectors = generate_tip_tail_vectors(transformed_field_array.copy(), sample_density_array, volume_box_array)

    generate_bild_file(tip_to_tail_vectors, transformed_field_array, args.ch, args.c, args.s, args.o)

main()
