import graphviz
from graphviz import digraph

dot = graphviz.Digraph('secondary RNA structure')
for char in SEQUENCE:
    dot.node(char, repr(char)) 