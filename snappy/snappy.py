def main():
    import argparse
    import sys
    import os
    from snappy.src.methods import parse_modkit, parse_fasta, add_contexts, construct_control_graph
    from snappy.src.main_stream import enrichment
    from snappy.src.generate_output import prepare_output, create_viz
    import numpy as np

    MIN_NUM_CONTEXTS = 30
    FRAC_MOD_THR = 80
    PERCENTILE = 0.1
    MIN_COV_THR = 10

    parser = argparse.ArgumentParser()

    parser.add_argument('-mk_bed', 
                        type=str, 
                        required=True, 
                        help='BED-file resulting from "modkit pileup"')
    parser.add_argument('-genome', 
                        type=str, 
                        required=True, 
                        help='Genome in FASTA format')
    parser.add_argument('-outdir',
                        type=str,
                        help='Output directory',
                        default='./Snappy_output')
    
    
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()


    out = args.outdir
    reference = parse_fasta(args.genome)

    try:        
        os.mkdir(out)
    
    except FileExistsError:

        print('The output directory already exists!')
        sys.exit()

    print('                       RUNNING ...                        ')
    print('===================================================================')

    modkit_df = parse_modkit(args.mk_bed)
    modkit_df = modkit_df.filter(modkit_df['n_valid'] > modkit_df['n_diff'])

    COV_THR = max(
        [int(np.percentile(modkit_df['n_valid'], PERCENTILE)), MIN_COV_THR]
    )
    print(f'Coverage threshold: {COV_THR}')

    modkit_df = modkit_df.filter(modkit_df['n_valid'] > COV_THR)
    modkit_df = add_contexts(modkit_df, reference)
    
    modkit_df = modkit_df.filter(modkit_df['context'].str.contains('.'*31))


    dict_freq_cntr, edge_count_cntr = construct_control_graph(reference)


    a6 = modkit_df.filter(modkit_df['mod_type'] == 'a')
    m5 = modkit_df.filter(modkit_df['mod_type'] == 'm')
    m4 = modkit_df.filter(modkit_df['mod_type'] == '21839')

    results = {}

    df = a6.filter(a6['frac_mod'] > FRAC_MOD_THR)
    results['a'] = enrichment(df, a6, dict_freq_cntr, 'A', MIN_NUM_CONTEXTS=MIN_NUM_CONTEXTS)


    df = m5.filter(m5['frac_mod'] > FRAC_MOD_THR)
    results['m'] = enrichment(df, m5, dict_freq_cntr, 'C', MIN_NUM_CONTEXTS=MIN_NUM_CONTEXTS)


    df = m4.filter(m4['frac_mod'] > FRAC_MOD_THR)
    results['21839'] = enrichment(df, m4, dict_freq_cntr, 'C', MIN_NUM_CONTEXTS=MIN_NUM_CONTEXTS)


    if len(results['a']) + len(results['m']) + len(results['21839']) == 0:
        print('Methylation motifs were not identified!')

    else:
        summary_table, results_table, regular_code = prepare_output(results, modkit_df)


        print('                       PREPARING THE OUTPUT                        ')
        print('===================================================================')
        
        summary_table.write_csv(f'{out}/Summary_table.tsv', separator="\t")
        results_table.write_csv(f'{out}/Results_table.tsv', separator="\t")

        # Create advanced out
        try:        
            os.mkdir(f'{out}/advanced_res')
        
        except FileExistsError:

            print('The output directory already exists!')

        regular_code.write_csv(f'{out}/advanced_res/Reg_exp.tsv', separator="\t")
        modkit_df.write_csv(f'{out}/advanced_res/Modkit_tab.tsv', separator="\t")

        os.mkdir(out + '/VIZ/')
        create_viz(results, modkit_df, out)

        
        print('                               DONE!                               ')
        print('===================================================================')

if __name__ == "__main__":
    main()