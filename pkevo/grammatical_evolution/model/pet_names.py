"""
This little script's only purpose is to generate silly names that can be attributed to 
the generated individuals. There is otherwise no real purpose to this script.
"""

import random

colors = [
    "red",
    "blue",
    "green",
    "yellow",
    "purple",
    "pink",
    "brown",
    "black",
    "white",
    "magenta",
    "violet"
]

fruits = [
    "apple", 
    "banana",
    "grape", 
    "strawberry", 
    "kiwi", 
    "mango", 
    "watermelon", 
    "blueberry", 
    "durian", 
    "lychee", 
    "jackfruit", 
    "pineapple",
    "mangosteen",
    "pomelo",
    "guava",
    "papaya",
    "nectarine",
    "raspberry"
]

fantasy_creatures = [
    "dragon",
    "phoenix",
    "griffin",
    "kelpie",
    "yeti",
    "banshee",
    "sphinx",
    "chimera",
    "wyvern",
    "hippogriff",
    "harpy",
    "cyclops",
    "basilisk",
    "gorgon",
    "kraken",
    "mermaid",
    "centaur",
    "pegasus",
    "siren",
    "striga",
    "djinn",
    "leprechaun",
    "unicorn",
    "gnome",
    "pixie"
]

def generate_random_name(creature_constraint=None):
    if creature_constraint != None:
        creature = creature_constraint
    else:
        creature = random.choice(fantasy_creatures)
    fruit = random.choice(fruits)
    color = random.choice(colors)
    name = f"{color}-{fruit}-{creature}"
    return name