from external_tools.mod_importer import mod

# Dictionary containing the mapping for rules to strategies
extended_strategy_mapping = dict()


# Define your custom (Sub-)strategy here
def user_custom_substrategy(strategy, pk_strategy_generator):
    # Define user's custom function logic here and alter the provided strategy
    print("I made it into my custom strategy function")
    return strategy


# Map your custom substrategy to a specific rule
#extended_strategy_mapping['<S>'] = user_custom_substrategy