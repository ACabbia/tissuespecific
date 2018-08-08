#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 17:18:33 2018

@author: acabbia
"""
import networkx as nx


def model_to_network(model, *args, **kwargs):
    """Convert a model into a networkx graph.
    Calls reactions_to_network with model.reactions.
    Parameters
    ----------
    model : cobra.Model
        The model.
    Returns
    -------
    networkx.MultiDiGraph
    See Also
    --------
    reactions_to_network
    """
    return reactions_to_network(model.reactions, *args, **kwargs)


def distance_based_on_molecular_formula(metabolite1, metabolite2, normalize=True):
    """Calculate the distance of two metabolites bases on the molecular formula
    Arguments
    ---------
    metabolite1 : Metabolite
        The first metabolite.
    metabolite2 : Metabolite
        The second metabolite.
    normalize : bool, optional
        If the distance should be normalized by the total number of elements in both metabolites (defaults to True).
    Returns
    -------
    float
        The distance between metabolite1 and metabolite2.
    """
    if len(metabolite1.elements) == 0 or len(metabolite2.elements) == 0:
        raise ValueError('Cannot calculate distance between metabolites %s and %s' % (metabolite1, metabolite2))
    elements = set(list(metabolite1.elements.keys()) + list(metabolite2.elements.keys()))
    distance = 0.
    for element in elements:
        distance += abs(metabolite1.elements.get(element, 0) - metabolite2.elements.get(element, 0))
    if normalize:
        return distance / sum(list(metabolite1.elements.values()) + list(metabolite2.elements.values()))
    else:
        return distance

def reactions_to_network(reactions, max_distance=0.3):
    """Convert a list of reactions into a networkx graph.
    Parameters
    ----------
    reactions : list
        The list of reactions.
    max_distance : float, optional
        A threshold on the normalized distance between two compounds. If distance is above this threshold,
        no edge between those compounds will be created.
    Returns
    -------
    networkx.MultiDiGraph
    See Also
    --------
    distance_based_on_molecular_formula
    """
    edges = list()
    for reaction in reactions:
        substrates = list(reaction.reactants)
        for substrate in substrates:
            products = list(reaction.products)
            for product in products:
                try:
                    distance = distance_based_on_molecular_formula(substrate, product, normalize=True)
                except ValueError:
                    distance = 0.
                if distance <= max_distance:
                    if reaction.reversibility:
                        edges.append((product, substrate, dict(reaction=reaction)))
                        edges.append((substrate, product, dict(reaction=reaction)))
                    elif reaction.lower_bound > 0:
                        edges.append((substrate, product, dict(reaction=reaction)))
                    else:
                        edges.append((product, substrate, dict(reaction=reaction)))

    multi_graph = nx.MultiDiGraph(edges)
    return multi_graph