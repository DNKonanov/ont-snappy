
from snappy.src.methods import get_re_from_degenerate
import polars as pl
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns 

def prepare_output(results, modkit_df):
    mod_decode = {'a' : '6mA', 'm' : '5mC', '21839' : '4mC'}
    summary_table = {"Motif": [], 
                    "Modification type": [], 
                    "Num positions": [], 
                    "Median FracMod": []}
    mods = []
    regular_code = {'Site' : [], 'RegExp' : []}
    
    for mod in results.keys():

        modkit_df_sub = modkit_df.filter(modkit_df['mod_type'] == mod)
        
        for site,_,_ in results[mod]:
            
            modkit_df_motif = modkit_df_sub.filter(modkit_df_sub['context'].str.contains(get_re_from_degenerate(site)))
            site_rename = f"{site[:15]}{mod_decode[mod]}{site[16:]}".strip(".*").replace('.', 'N')
            
            regular_code['Site'].append(site_rename)
            regular_code['RegExp'].append(get_re_from_degenerate(site))
            
            summary_table["Motif"].append(site_rename)
            summary_table["Modification type"].append(mod_decode[mod])
            summary_table["Num positions"].append(len(modkit_df_motif))
            summary_table["Median FracMod"].append(np.round(np.median(list(modkit_df_motif['frac_mod'])), 2))
            modkit_df_motif = modkit_df_motif.with_columns(
                                        pl.Series(name="motif", values=[site_rename for i in range(len(modkit_df_motif))]),
                                        )
            modkit_df_motif = modkit_df_motif.with_columns(
                                        pl.Series(name="new_mod", values=[mod_decode[mod] for i in range(len(modkit_df_motif))]),
                                        )
            
            mods.append(modkit_df_motif)
    # Summary table
    summary_table = pl.from_dict(summary_table)
    mods = pl.concat(mods, how="align")

    new_names = {'contig' : 'Contig', 
                'start' : 'Coordinate', 
                'strand' : 'Strand', 
                'new_mod' : 'Modification', 
                'n_valid' : 'Coverage', 
                'motif' : 'Motif'}

    results_table =  mods[['contig', 'start', 'strand', 'new_mod', 'frac_mod', 'n_valid', 'motif']].rename(new_names)
    
    regular_code = pl.from_dicts(regular_code)

    return summary_table, results_table, regular_code

def create_viz(results, modkit_df, out):
    
    mod_decode = {'a' : '6mA', 'm' : '5mC', '21839' : '4mC'}

    for mod in results.keys():

        modkit_df_sub = modkit_df.filter(modkit_df['mod_type'] == mod)
        
        for site,_,_ in results[mod]:
        
            modkit_df_motif = modkit_df_sub.filter(modkit_df_sub['context'].str.contains(get_re_from_degenerate(site)))
            site_rename = f"{site[:15]}{mod_decode[mod]}{site[16:]}".strip(".*").replace('.', 'N')
            num_contig = len(set(modkit_df_motif['contig']))

            fig, axs = plt.subplots(num_contig, 1, figsize=(12, len(set(modkit_df_motif['contig']))*2))
            c = 0
            
            for contig in np.sort(list(set(modkit_df_motif['contig']))):
                if num_contig == 1: 
                    ax=axs
                else:
                    ax=axs[c]

                contig_subset = modkit_df_motif.filter(modkit_df_motif['contig'] == contig)
                sns.lineplot(x=contig_subset['start'], y=contig_subset['frac_mod'], linewidth=.5, c='grey', ax=ax)
                ax.set_ylabel('Fraction modified, %')
                ax.set_xlabel(f'Position in {contig}')
                ax.set_title(contig)
                ax.set_ylim(0, 100)
                ax.set_xlim(np.min(list(contig_subset['start']))-100, np.max(list(contig_subset['start']))+100)
                ax.spines['right'].set_visible(False)
                ax.spines['top'].set_visible(False)
                c += 1
            fig.suptitle(site_rename, fontsize=14)
            plt.tight_layout()
            plt.savefig(f'{out}/VIZ/{site_rename}_localization.pdf', bbox_inches='tight')

            fig, axs = plt.subplots(1, 1)
            sns.histplot(modkit_df_motif['frac_mod'], ax=axs, color='grey')
            axs.set_title(site_rename, fontsize=14)
            axs.set_xlabel('Fraction modified, %')
            axs.spines['right'].set_visible(False)
            axs.spines['top'].set_visible(False)
            sns.despine(trim=True)
            plt.savefig(f'{out}/VIZ/{site_rename}_FracMod.pdf', bbox_inches='tight')