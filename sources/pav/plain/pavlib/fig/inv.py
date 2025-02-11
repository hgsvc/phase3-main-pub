"""
Inversion figures
"""

import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt

import pavlib
import kanapy


def kde_density_base(df_kde, region_qry, width=7, height=4, dpi=300, flank_whiskers=False):
    """
    Base k-mer density plot using a k-mer density DataFrame.

    :param df_kde: Inversion call DataFrame.
    :param region_qry: Query region plot is generated over. Must match the region in `df_kde`.
    :param width: Figure width.
    :param height: Figure height.
    :param dpi: Figure DPI.
    :param flank_whiskers: If `True`, show whiskers above or below points indicating if they match the upstream or
        downstream flanking inverted duplication.

    :return: Plot figure object.
    """

    # Make figure
    fig = plt.figure(figsize=(width, height), dpi=dpi)

    ax1, ax2 = fig.subplots(2, 1)


    ## Smoothed state (top pane) ##

    # Flanking DUP whiskers
    if flank_whiskers:
        raise NotImplementedError('Flanking whiskers are not yet implemented.')
        # for index, subdf in df_kde.loc[
        #     ~ pd.isnull(df_kde['MATCH'])
        # ].groupby(
        #     ['STATE', 'MATCH']
        # ):
        #
        #     # Whisker y-values
        #     ymin, ymax = sorted(
        #         [
        #             -1 * (subdf.iloc[0]['STATE_MER'] - 1),
        #             -1 * (subdf.iloc[0]['STATE_MER'] - 1) + (0.25 if index[1] == 'SAME' else -0.25)
        #         ]
        #     )
        #
        #     # Plot whiskers
        #     ax1.vlines(
        #         x=subdf['INDEX'] + region_qry.pos,
        #         ymin=ymin, ymax=ymax,
        #         color='dimgray',
        #         linestyles='solid',
        #         linewidth=0.5
        #     )

    # Points (top pane)
    ax1.scatter(
        x=np.asarray(df_kde.loc[df_kde['STATE_MER'] == 0, 'INDEX']) + region_qry.pos,
        y=np.repeat(1, np.sum(df_kde['STATE_MER'] == 0)),
        color='blue',
        alpha=0.2
    )

    ax1.scatter(
        x=np.asarray(df_kde.loc[df_kde['STATE_MER'] == 1, 'INDEX']) + region_qry.pos,
        y=np.repeat(0, np.sum(df_kde['STATE_MER'] == 1)),
        color='purple',
        alpha=0.2
    )

    ax1.scatter(
        x=np.asarray(df_kde.loc[df_kde['STATE_MER'] == 2, 'INDEX']) + region_qry.pos,
        y=np.repeat(-1, np.sum(df_kde['STATE_MER'] == 2)),
        color='red',
        alpha=0.2
    )

    # Max density line (smoothed state call)
    ax1.plot(
        df_kde['INDEX'] + region_qry.pos,
        (df_kde['STATE'] - 1) * -1,
        color='black'
    )

    # Plot aestetics
    ax1.get_xaxis().set_major_formatter(
        mpl.ticker.FuncFormatter(lambda x, p: format(int(x), ','))
    )

    ax1.set_yticks(np.asarray([-1, 0, 1]))
    ax1.set_yticklabels(np.array(['Rev', 'Fwd+Rev', 'Fwd']))


    ## Density (bottom pane) ##

    ax2.plot(
        df_kde['INDEX'] + region_qry.pos,
        df_kde['KDE_FWD'],
        color='blue'
    )

    ax2.plot(
        df_kde['INDEX'] + region_qry.pos,
        df_kde['KDE_FWDREV'],
        color='purple'
    )

    ax2.plot(
        df_kde['INDEX'] + region_qry.pos,
        df_kde['KDE_REV'],
        color='red'
    )

    # Plot aestetics

    ax2.get_xaxis().set_major_formatter(
        mpl.ticker.FuncFormatter(lambda x, p: format(int(x), ','))
    )

    ax2.set_xlabel('{} ({:,d} - {:,d})'.format(
        region_qry.chrom,
        region_qry.pos + 1,
        region_qry.end
    ))

    ax1.tick_params(labelbottom=False)

    for label in ax2.get_xticklabels():
        label.set_rotation(30)
        label.set_ha('right')


    ## Return plot and axes ##
    return fig


def dotplot_inv_call(
    inv_call, ref_fa, qry_fa=None, seq_qry=None,
    region_ref=None, region_qry=None,
    k=32
):
    """
    Make a dotplot of an inversion call.

    :param inv_call: pavlib.inv.InvCall object describing the inversion.
    :param ref_fa: Reference FASTA.
    :param qry_fa: Query FASTA file name.
    :param seq_qry: query alignment sequence or `None`. If `None`, then `qry_fa` must be set.
    :param region_ref: Reference region. If None, uses inv_call.region_ref_discovery.
    :param region_qry: Query region. If none, uses inv_call.region_tig_discovery.
    :param k: K-mer size.

    :return: A plot object. Before it is discarded, this object should be closed with `matplotlib.pyplot.close()` to
        free memory.
    """

    if region_ref is None:
        region_ref = inv_call.region_ref

    if region_qry is None:
        region_qry = inv_call.region_qry

    # Get reference sequence
    seq_ref = pavlib.seq.region_seq_fasta(region_ref, ref_fa, False)

    # Get contig sequence
    if seq_qry is None:

        if qry_fa is None:
            raise RuntimeError('Cannot get contig sequence: tig_fa is None')

        seq_qry = pavlib.seq.region_seq_fasta(region_qry, qry_fa)

    # Create plot config
    plot_config = {
        'label_x': '{} ({:,d} - {:,d})'.format(region_qry.chrom, region_qry.pos + 1, region_qry.end),
        'label_y': '{} ({:,d} - {:,d})'.format(region_ref.chrom, region_ref.pos + 1, region_ref.end),
        'start_x': region_qry.pos,
        'start_y': region_ref.pos,
        'k': k
    }

    # Add annotations
    anno_list = list()

    if inv_call.region_ref_outer is not None and inv_call.region_ref_outer != inv_call.region_ref:
        anno_list.append(
            {
                'type': 'hshade',
                'y1': inv_call.region_ref_outer.pos + 1,  # From BED to 1-based closed coordinates
                'y2': inv_call.region_ref_outer.end,
                'color': 'lightseagreen'
            }
        )

    if inv_call.region_qry_outer is not None and inv_call.region_qry_outer != inv_call.region_qry:
        anno_list.append(
            {
                'type': 'vline',
                'x': np.array([inv_call.region_qry_outer.pos, inv_call.region_qry_outer.end]),
                'ymin': region_ref.pos + 1,
                'ymax': region_ref.end,
                'color': 'lightseagreen',
                'style': 'solid',
                'alpha': 0.4
            }
        )

        anno_list.append(
            {
                'type': 'hshade',
                'y1': inv_call.region_ref.pos + 1,  # From BED to 1-based closed coordinates
                'y2': inv_call.region_ref.end,
                'color': 'seagreen'
            }
        )

        anno_list.append(
            {
                'type': 'vline',
                'x': np.array([inv_call.region_qry.pos, inv_call.region_qry.end]),
                'ymin': region_ref.pos + 1,
                'ymax': region_ref.end,
                'color': 'seagreen',
                'style': 'dashed',
                'alpha': 0.4
            }
        )

    # Make plot
    fig = kanapy.plot.dotplot.dotplot(
        seq_x=seq_qry,
        seq_y=seq_ref,
        config=plot_config,
        title=inv_call.id,
        anno_list=anno_list
    )

    # Return figure
    return fig
