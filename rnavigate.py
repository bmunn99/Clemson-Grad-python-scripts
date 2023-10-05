# https://rnavigate.readthedocs.io/en/latest/
# Import rnavigate
import rnavigate as rnav
import numpy as np

# Load in your RNA data (ct and shape file - must be profile.txt file)
AGO3_WT = rnav.Sample(
    sample="AGO3_WT",
    fasta='/data2/lackey_lab/shalam/rnastructure/shapemapper/FASTA/full/AGO3_WT.fa',
    ct='/data2/lackey_lab/shalam/rnastructure/fold/CT/Nextseq/exp/3_WT.ct',
    shapemap='/data2/lackey_lab/shalam/rnastructure/shapemapper/AGO3_results/Nextseq/3_WT/AGO3_WT_AGO3_WT_CDS_profile.txt',
    pairprob='/data2/lackey_lab/shalam/rnastructure/partition/dp/exp/AGO3_WT.dp',
    annotations={
        'seq_source': 'fasta',
        'E638A_mut': {
            'spans': [[1912, 1914]],
            'color': 'green',
        },
        'P2_mut': {
            'spans': [[673, 675]],
            'color': 'green',
        },
        'P1_mut': {
            'spans': [[1519, 1521]],
            'color': 'green',
        }
    })

AGO3_E638A = rnav.Sample(
    sample="AGO3_E638A",
    fasta='/data2/lackey_lab/shalam/rnastructure/shapemapper/FASTA/full/AGO3_E638A.fa',
    ct='/data2/lackey_lab/shalam/rnastructure/fold/CT/Nextseq/exp/E638A.ct',
    shapemap='/data2/lackey_lab/shalam/rnastructure/shapemapper/AGO3_results/Nextseq/E638A/AGO3_E638A_AGO3_E638A_CDS_profile.txt',
    pairprob='/data2/lackey_lab/shalam/rnastructure/partition/dp/exp/AGO3_E638A.dp',
    annotations={
        'seq_source': 'fasta',
        'E638A_mut': {
            'spans': [[1912, 1914]],
            'color': 'green',
        }
    })

AGO3_P2 = rnav.Sample(
    sample="AGO3_P2",
    fasta='/data2/lackey_lab/shalam/rnastructure/shapemapper/FASTA/full/AGO3_P2.fa',
    ct='/data2/lackey_lab/shalam/rnastructure/fold/CT/Nextseq/exp/P2.ct',
    shapemap='/data2/lackey_lab/shalam/rnastructure/shapemapper/AGO3_results/Nextseq/P2/AGO3_P2_AGO3_P2_CDS_profile.txt',
    pairprob='/data2/lackey_lab/shalam/rnastructure/partition/dp/exp/AGO3_P2.dp',
    annotations={
        'seq_source': 'fasta',
        'P2_mut': {
            'spans': [[673, 675]],
            'color': 'green',
        }
    })

AGO3_P1 = rnav.Sample(
    sample="AGO3_P1",
    ct='/data2/lackey_lab/shalam/rnastructure/fold/CT/Nextseq/exp/P1.ct',
    fasta='/data2/lackey_lab/shalam/rnastructure/shapemapper/FASTA/full/AGO3_P1.fa',
    shapemap='/data2/lackey_lab/shalam/rnastructure/shapemapper/AGO3_results/Nextseq/P1/AGO3_P1_AGO3_P1_CDS_profile.txt',
    pairprob='/data2/lackey_lab/shalam/rnastructure/partition/dp/exp/AGO3_P1.dp',
    annotations={
        'seq_source': 'fasta',
        'P1_mut': {
            'spans': [[1519, 1521]],
            'color': 'green',
        }
    })

AGO2_WT = rnav.Sample(
    sample="AGO2_WT",
    ct='/data2/lackey_lab/shalam/rnastructure/fold/CT/Nextseq/exp/2_WT.ct',
    fasta='/data2/lackey_lab/shalam/rnastructure/shapemapper/FASTA/full/AGO2_WT.fa',
    shapemap='/data2/lackey_lab/shalam/rnastructure/shapemapper/AGO3_results/Nextseq/2_WT/AGO2_WT_AGO2_WT_CDS_profile.txt',
    pairprob='/data2/lackey_lab/shalam/rnastructure/partition/dp/exp/AGO2_WT.dp')

# Make an arc diagram comparison
WT_E638A = rnav.plot_arcs_compare(samples=[AGO3_WT, AGO3_E638A], ct=None, seq_source='fasta', 
                                 interactions='pairprob',
                                  interactions_filter={'Probability_ge': 0.4})
WT_P2 = rnav.plot_arcs_compare(samples=[AGO3_WT, AGO3_P2], ct=None, seq_source='fasta', 
                                 interactions='pairprob',
                                  interactions_filter={'Probability_ge': 0.4})
WT_P1 = rnav.plot_arcs_compare(samples=[AGO3_WT, AGO3_P1], ct=None, seq_source='fasta', 
                                 interactions='pairprob',
                                  interactions_filter={'Probability_ge': 0.4})

WT_E638A.save('WT_E638A.png')
WT_P2.save('WT_P2.png')
WT_P1.save('WT_P1.png')

# Make an arc diagram comparison zoomed into region of interest
WT_E638A = rnav.plot_arcs(samples=[AGO3_WT], region=[1662, 2162], annotations=['E638A_mut'], ct=None,
                          seq_source='fasta')
rnav.fit_data_list(AGO3_E638A, ['ct'], AGO3_WT.data['ct'])
WT_E638A.add_patches(
    ax=WT_E638A.axes[0, 0],
    data=AGO3_E638A.data['pairprob'],
    panel='bottom',
    annotation_gap=4)
WT_E638A.add_patches(
    ax=WT_E638A.axes[0, 0],
    data=AGO3_WT.data['pairprob'],
    panel='top',
    annotation_gap=4)

WT_E638A.save('WT_E638A_zoom.png')
