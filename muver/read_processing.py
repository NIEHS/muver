def remove_diff_chr_pairs(in_sam, out_sam):
    '''
    Go through a SAM file and remove entries where the other read is on a
    different chromosome.
    '''
    with open(in_sam) as in_file, open(out_sam, 'w') as OUTPUT:

        for line in in_file:

            if line[0] == '@':
                OUTPUT.write(line)

            else:
                line_split = line.strip().split('\t')
                if line_split[6] == '=':
                    OUTPUT.write(line)
