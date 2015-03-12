import re

class graph(object):
    def __init__(self, V = None, E = None):
        if V == None:
            self.verts = set()
        elif type(V) in set([list, set, tuple]):
            self.verts = set(V)
        else:
            raise ValueError, 'bad V value'
            
        if E == None:
            self.edges = set()
        elif type(V) in set([list, set, tuple]):
            self.edges = set(E)
        else:
            raise ValueError, 'bad E value'
    
    def __iter__(self):
        for e in self.edges:
            yield e
            
    def __repr__(self):
        vertStr = ','.join(str(v) for v in self.verts)
        edgeStr = ','.join(str(e) for e in self)
        return 'V = ' + vertStr + '\nE = ' + edgeStr
        
    def __len__(self):
        return len(self.edges)
        
    def add(self, e):
        if type(e) == tuple:
            self.edges.add(e)
            
    def addVert(self, v):
        self.verts.add(v)
        
class weightedGraph(graph):
    def __init__(self, G = None, W = None):
        if G == None:
            self.verts = set()
            self.edges = dict()
        elif type(G) == graph:
            self.verts = set(G.verts)
            self.edges = {e:1 for e in G}
        else:
            raise ValueError, "G must be a graph in initialization"
        
        if W != None:
            for e, w in W:
                self.edges[e] = w
                
    def add(self, e, w):
        self.edges[e] = w
                
    def __iter__(self):
        return self.edges.iteritems()
                
    def __repr__(self):
        vertStr = ','.join(str(v) for v in self.verts)
        edgeStr = ','.join(str(e) +':'+ str(w) for e,w in self)
        return 'V = ' + vertStr + '\nE = ' + edgeStr
        
           
if __name__ == "__main__":
    G = graph(V = set([1,2]))
    G.edges.add((1,2))
    G.edges.add((1,3))
    wG = weightedGraph(G)
    