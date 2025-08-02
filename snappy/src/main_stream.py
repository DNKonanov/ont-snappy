
from snappy.src.methods import get_dict_freq_exp, get_chi_stats, get_initial_max_edge, get_re_template, get_best_positions, curate_motif, get_re_from_degenerate
from scipy.stats import chi2_contingency

def get_chi_stat_ind(template, df, df_cntr):

    sdf = df.filter(df['context'].str.contains(get_re_from_degenerate(template)))
    sdf_cntr = df_cntr.filter(df_cntr['context'].str.contains(get_re_from_degenerate(template)))

    f1, s1, f2, s2 = len(sdf), len(df), len(sdf_cntr), len(df_cntr)

    return chi2_contingency(
        [[f1, s1], [f2, s2]]
    )



def enrichment(df, df_cntr, dict_freq_cntr, central_base, MIN_NUM_CONTEXTS=50, OFFSETS_RANGE=list(range(12)), MAX_MOTIFS=100):

    found_templates = []
    found_templates_cnt = []
    
    df_ = df.clone()

    iteration = 1
    print('\n')
    print(f'-----ITERATION {iteration}: {len(df)} unexplained contexts in the data-----')


    for i in range(MAX_MOTIFS):
        if len(df) < MIN_NUM_CONTEXTS:
            break
        dict_freq_exp = get_dict_freq_exp(df, OFFSETS_RANGE)
        chi_stat_tables = get_chi_stats(dict_freq_exp, dict_freq_cntr)
        initial_max_edge = get_initial_max_edge(chi_stat_tables)
        print('initial max edge:', initial_max_edge)

        re_template = get_re_template(*initial_max_edge[0], initial_max_edge[2])

        motif_positions = get_best_positions(re_template, df['context'])
        extended_re_template = '.'*(motif_positions[0]) + re_template + '.'*(31 - motif_positions[1])
        
        sdf = df.filter(df['context'].str.contains(re_template))
        if len(sdf) < MIN_NUM_CONTEXTS: 
            df = df.filter(df['context'].str.contains(get_re_from_degenerate(re_template)[1:-1]) == False)            
            continue

        sig = False
        for t in found_templates:
            if t.replace('.', '') in extended_re_template:
                print('Skipped as a submotif!')
                sig = True
                break
        
        if sig == False:
            curated_template = curate_motif(extended_re_template, central_base, df_cntr)

            
            sdf = df.filter(df['context'].str.contains(get_re_from_degenerate(curated_template)))
            sdf_cntr = df_cntr.filter(df_cntr['context'].str.contains(get_re_from_degenerate(curated_template)))

            f1, s1, f2, s2 = len(sdf), len(df), len(sdf_cntr), len(df_cntr)

            matrix = [[f1, s1], [f2, s2]]
            
            
            
            chi_res = chi2_contingency(
                matrix
            )

            print('Confidence:', round(chi_res[0], 1))
            
            if matrix[0][0] > MIN_NUM_CONTEXTS and chi_res[0] > 100:    
                found_templates.append(curated_template.replace('N', '.'))
                
                found_templates_cnt.append(
                    (
                        curated_template.replace('N', '.'), 
                        f1,
                        chi_res[0],
                    )
                )
                
            

        

        
        prev_len = len(df)
        df = df.filter(df['context'].str.contains(get_re_from_degenerate(curated_template)[1:-1]) == False)
        post_len = len(df)
        if prev_len == post_len:
            df = df.filter(df['context'].str.contains(get_re_from_degenerate(re_template)[1:-1]) == False)
        

        if len(df) < MIN_NUM_CONTEXTS:
            break
        print('\n')
        iteration += 1
        print(f'-----ITERATION {iteration}: {len(df)} unexplained contexts in the data-----')

    return found_templates_cnt