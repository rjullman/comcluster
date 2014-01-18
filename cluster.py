from __future__ import print_function
from __future__ import division
from graph_tool.all import *
from argparse import ArgumentParser
from unidecode import unidecode
from random import random
import simplejson
import itertools
import sys

# ------------------------------------------------------------------------------

def unidecode_safe(ustring):
    """Returns a string from a unicode string failing gracefully in corner cases."""
    return unidecode(ustring) if isinstance(ustring,unicode) else ustring

def add_edge_props(g, e, title, year, citation):
    """Sets the edge properties for an edge.

    Keyword Arguments:
    g        -- graph containing the edge to set the property
    e        -- the edge for which to set the properties
    title    -- title of the paper on that edge
    year     -- year of publication/submission
    citation -- boolean indication if edge contains a citation

    """
    title = unidecode_safe(title)
    year = unidecode_safe(year) 

    g.edge_properties["title"][e] = title
    g.edge_properties["year"][e] = year
    g.edge_properties["citation"][e] = citation

def get_vertex_by_author(g, author):
    """Returns the vertex in graph g corresponding to the given author.

    If a vertex with that author does not exist a new one is created,
    added to the graph and returned.

    Keyword Argument:
    g      -- graph in which to search for the author
    author -- the string or unicode string of the author's name

    """
    author = unidecode_safe(author)
    author_to_vertex = g.graph_properties['authorToVertex']

    if author in author_to_vertex:
        v = g.vertex(author_to_vertex[author])
    else:

        global alert
        if alert:
            print(author)
        
        v = g.add_vertex()
        g.vertex_properties['author'][v] = author
        author_to_vertex[author] = int(v)

    return v
        
def add_paper_to_graph(g, paper_data, citation_weight):
    """Adds the paper data to graph g.

    Keyword Arguments:
    g                -- graph to add the data to
    paper_data       -- dictionary containing the paper data
    citattion_weight -- weight to give to a citation edge

    The paper data must contain at least a list of authors in order to
    be added to the graph.
   
    Lower citation weight means that the connected authors are more
    related.

    """

    # Not enough data to add the paper to the graph
    if not 'author' in paper_data:
        return

    # Extract paper data
    authors = paper_data['author'] if isinstance(paper_data['author'], list) else [paper_data['author']]
    title = paper_data['title'] if 'title' in paper_data else None
    year = paper_data['year'] if 'year' in paper_data else None
    citations = paper_data['citations'] if 'citations' in paper_data else None
    
    # Edge weight property
    edge_dists = g.edge_properties['dists']

    # Add paper with single author
    if len(authors) is 1:
        v = get_vertex_by_author(g, authors[0])
        e = g.add_edge(v,v)
        add_edge_props(g, e, title, year, False)
        edge_dists[e] = 1
        g.vertex_properties["num_papers"][v] += 1

    # Add paper with co-authors
    else:
        for i, a1 in enumerate(authors):
            v1 = get_vertex_by_author(g, a1)
            g.vertex_properties["num_papers"][v1] += 1
            for a2 in authors[i+1:]:
                v2 = get_vertex_by_author(g, a2)
                e = g.add_edge(v1,v2)
                add_edge_props(g, e, title, year, False)
                edge_dists[e] = 60.0 / len(authors)

    # Add citations
    if citations:
        for a in authors:
            for c in citations:
                v1 = get_vertex_by_author(g, a)
                v2 = get_vertex_by_author(g, c)
                e = g.add_edge(v1,v2)
                add_edge_props(g, e, None, None, True)
                edge_dists[e] = citation_weight

# ------------------------------------------------------------------------------

def print_partition(g, partition):
    """Displays the given partition of graph g to stdout.
    
    Each vertex array in the partition is sorted in descending order
    of influence in the graph.

    """
    for vertices in partition:
        sort_vertices(g,vertices)
        np = map(lambda v : g.vertex_properties['author'][v], vertices)
        print(np)
        print()

def vertices_to_names(g, vertices):
    return map(lambda v : g.vertex_properties['author'][v], vertices)

def degree(g, v):
    """Returns the degree of v in g"""
    return g.vertex_properties["num_papers"][v]

def sort_vertices(g, vertices):
    """Orders the vertex set in graph g"""
    return list.sort(vertices, key = lambda v : degree(g,v), reverse=True)

# ------------------------------------------------------------------------------

def subgraph_filter(g, vertices):
    """Filters graph g to contain only selected vertices."""
    fprop = g.new_vertex_property("bool")
    for v in vertices:
        fprop[v] = True
    g.set_vertex_filter(None)
    g.set_vertex_filter(fprop)
    return g

def partition_from_map(partition_map):
    """Returns a partition (array of vertex arrays) from a graph_tool property map."""
    partition = [[] for i in range(g.num_vertices())]
    for v in g.vertices():
        partition[partition_map[v]].append(v)
    return filter(lambda x : any(x), partition)

def gen_partition(g, size_cutoff = -1, depth = 1, verbose=False):
    """Returns an array of vertice arrays that partition graph g.

    Keyword arguments:
    g          -- the graph to partition
    size_cutff -- the cluster size below which to stop subdividing clusters (default -1, no cutoff)
    depth      -- the iteration depth past which to stop subdividing clusters (default 1)

    """
    if verbose:
        print("Processing cluster ({0} vertices - depth {1})"
              .format(g.num_vertices(), depth), file=sys.stderr)

    if g.num_vertices() < size_cutoff:
        return [list(g.vertices())]

    # Generate partitions of the given graph
    partition_map = minimize_blockmodel_dl(g, eweight = g.edge_properties['dists'])[0]
    partition = partition_from_map(partition_map)

    if depth > 1 :
        # Recurse by re-partitioning the graph restricted to each prior partition
        partition = list(itertools.chain(
            *map(lambda vs : 
                 gen_partition(subgraph_filter(g,vs), size_cutoff, depth - 1, verbose),
                 partition_from_map(partition_map))))

    return partition

# ------------------------------------------------------------------------------

def reps_per_partition(tot_reps, weights):
    """Returns a integral partitions (integer list) of tot_reps by weights

    Note that the sum of the returned may be slightly larger (but not
    smaller) than tot_reps.

    tot_reps -- the number of representatives to distrubutes
    weights  -- a list of percentages of deserved by each partition

    """
    # Representatives per partition
    repspp = [0 for i in range(len(weights))]
    deserving = [0 for i in range(len(weights))]

    assigned_reps = 0
    for i in range(len(repspp)):
        ar = int(round(weights[i] * tot_reps))
        assigned_reps += ar
        repspp[i] = ar
        # Compute how much of an extra rep each partition deserves
        deserving[i] = weights[i] * tot_reps - ar

    if assigned_reps >= tot_reps:
        return repspp

    most_deserving = sorted(range(len(deserving)), key=lambda k: deserving[k], reverse=True)
    for i in most_deserving[0: tot_reps - assigned_reps]:
        assert (deserving[i] < 1), "Partition deserves more than a fractional representative?"
        assigned_reps += 1
        repspp[i] += 1

    assert (assigned_reps >= tot_reps), "Did not assign all representatives?"
    return repspp

def same_name(a_sim, a_com):
    a_sim = a_sim.lower();
    a_com = a_com.lower();
    for part in a_sim.split(' '):
        if not part in a_com:
            return False
    return True

def author_list_contains(authors, author):
    for a in authors:
        if same_name(a, author):
            return True
    return False

def get_reps(partition, tot_reps, num_alt, query_file, filter_file):
    # Total number of authors
    tweight = reduce(lambda a, l : a + len(l), partition, 0)
    # Start by weighting partitions by how many members they have
    pweight = map(lambda l : len(l)/tweight, partition)

    if query_file:
        authors = map(lambda s : s.strip(), open(args.query_file, 'r').read().split(','))
        ### do a double loop to collect found authors
        # Only consider partitions with authors in the query file
        partition = filter(lambda l : any(set(l) & set(authors)), partition)
        # Weight each partition by the number of querried authors they contain
        pweight = map(lambda l : len(set(l) & set(authors))/len(authors), partition)

    if filter_file:
        # Remove the desired authors from the representatives pool
        authors = map(lambda s : s.strip(), open(args.filter_file, 'r').read().split(','))
        partition = map(lambda l : [x for x in l if not author_list_contains(authors,x)], partition)

    # Return a list of representatives for the partition
    reps_per_part = reps_per_partition(tot_reps, pweight)
    reps = itertools.imap(lambda l, r : l[0:r + num_alt], partition, reps_per_part)
    reps = filter(lambda l : len(l) > num_alt, reps)
    return reps

# ------------------------------------------------------------------------------

def draw(g, reps, partition, draw_file, max_iters):
    # Remove all proir filtering of the graph
    g.set_vertex_filter(None)
    block = g.new_vertex_property("vector<double>")
    colors = []

    # Create array of possible colors
    for v in g.vertices():
        colors.append([random(), random(), random(), 0.9])

    # Set colors for each vertex depending on block
    for v in g.vertices():
        for i in range(len(partition)):
            if v in partition[i]:
                block[v] = colors[i]
                break

    # Make the representative vertices slightly larger
    size = g.new_vertex_property("int")
    for v in g.vertices():
        size[v] = 5
    for rep in reps:
        v = get_vertex_by_author(g, rep)
        size[v] = 8

    # Compute the arf layout of the graph and output it to the specified file
    pos = arf_layout(g, max_iter=max_iters)
    graph_draw(g, pos=pos, vertex_size=size, vertex_fill_color=block, output=draw_file)

def rep_dist(g, reps):
    g.set_vertex_filter(None)
    avg_dist = 0
    disconnected = 0

    for src in g.vertices():
        dists = shortest_distance(g, source=src)
    
        min_dist = sys.maxint
        for rep in reps:
            v = get_vertex_by_author(g, rep)
            if dists[v] < min_dist:
                min_dist = dists[v]

        if min_dist > g.num_edges():
            disconnected = disconnected + 1
            min_dist = 0

        avg_dist += min_dist

    num_vert = g.num_vertices() - disconnected - len(reps)
    if num_vert <= 0:
        avg_dist = 0
    else:
        avg_dist /= num_vert

    return avg_dist, disconnected

# ------------------------------------------------------------------------------

# Setup command line parser
p = ArgumentParser(description="Clusters authors into communities by coauthorship graph.")
p.add_argument('files', nargs='+', help='the files containing bibliographical data in json format')
p.add_argument('-n', dest='tot_reps', type=int, default=8, help='the number of people to select to represent the network, default 8')
p.add_argument('-i', dest='iter', type=int, default=3, help='the number of times to recursively sub-divide clusters, default 3')
p.add_argument('-c', dest='cutoff', type=int, default=50, help='the cluster size to stop further cluster sub-division, default 50')
p.add_argument('-w', dest='cweight', type=float, default=1, help='the weight of a citation edge (lower means more connected), default 1')
p.add_argument('-q', dest='query_file', help='the file containing a (comma delimited) list of authors for which to find representatives for (restricts graph data just to nodes around these authors).')
p.add_argument('-f', dest='filter_file', help='the file containing a (comma delimited) list of authors to not select as representatives of the graph')
p.add_argument('-a', dest='num_alt', type=int, default=0, help='the number of alternate representatives to choose for each group')
p.add_argument('-o', dest='draw_file', help='when given a graphical representation of the graph will be output to this file, see graph-tool graph_draw() for supported formats')
p.add_argument('-di', dest='draw_iters', type=int, default=250, help='when drawing output filename is given this determines the max number of iterations for arf layout (see graph-tool arf_layout) of the graph, default 250. When set to 0 the algorithm will iterate until a stable state is reached.')
p.add_argument('--no-stats', dest='stats', default=True, action='store_false', help='do not display stats that evaluate the quality of the selected representatives for the community')
p.add_argument('-s', dest='author_file', help='using the authors (comma delimited) as the representatives computes and displays stats that evaluates the quality that they represent the community')
p.add_argument('-v', dest='verbose', action='store_true', help='show verbose logging')
args = p.parse_args()

# Load in data from file
json_data_list = itertools.chain(
    *map(lambda f :
         map(lambda x : simplejson.loads(x), 
             open(f, "r").read().strip().split('\n')),
         args.files))

# Define meta data on the graph
g = Graph(directed=False)
g.edge_properties["title"] = g.new_edge_property("string")
g.edge_properties["year"] = g.new_edge_property("string")
g.edge_properties["citation"] = g.new_edge_property("bool")
g.edge_properties["dists"] = g.new_edge_property("double")
g.vertex_properties["author"] = g.new_vertex_property("string")
g.vertex_properties["num_papers"] = g.new_vertex_property("double")
g.graph_properties["authorToVertex"] = g.new_graph_property("object")
g.graph_properties["authorToVertex"] = {};

alert = False

# Create the graph
map(lambda x : add_paper_to_graph(g, x, args.cweight), json_data_list)

alert = True

# Read the community representatives from a file
if args.author_file:
    reps = map(lambda s : unidecode_safe(s.strip()), open(args.author_file, 'r').read().split(','))
    for rep in reps:
        print(rep)

# Otherwise compute the representatives for the community
else:
    # Partition the graph
    partition = gen_partition(g, size_cutoff=args.cutoff, depth=args.iter, verbose=args.verbose)

    # Sort the partition and convert the vertex lists to names of authors
    for l in partition:
        sort_vertices(g,l)
    author_partition = map(lambda l : vertices_to_names(g,l), partition)

    # Output the reps
    rep_lists = get_reps(author_partition, args.tot_reps, args.num_alt, args.query_file, args.filter_file)
    print("\n--------------------\n")
    for reps in rep_lists:
        for idx, author in enumerate(reps):
            if idx < len(reps) - args.num_alt:
                print(author)
            else:
                print("A:", author)
        print("\n--------------------\n")

    # flatten the representatives list
    reps = list(itertools.chain(*rep_lists))

# Display quality of selected representatives
if args.stats or args.author_file:
    avg_dist, disconnected = rep_dist(g, reps)
    print("Vertex Count:", g.num_vertices())
    print("Edge Count:", g.num_edges())
    print("Number of Reps:", len(reps))
    print("Avg Vertex Rep Distance:", avg_dist)    
    print("Unrepresented Vertices:", disconnected)

# Output a visual representation
if args.draw_file:
    draw(g, reps, partition, args.draw_file, args.draw_iters)
