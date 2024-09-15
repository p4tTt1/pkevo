import logging
import os

# from https://stackoverflow.com/questions/50981906/change-default-location-log-file-generated-by-logger-in-python#answer-69693313
# and https://docs.python.org/3/howto/logging-cookbook.html#using-logging-in-multiple-modules
def get_logger(module_name, log_file = "logs.log", log_level = logging.DEBUG):
    

    current_dir = os.path.dirname(os.path.abspath(__file__)) # should be /pkevo/config
    log_filepath = os.path.join(current_dir, "../..", "logs", log_file) # should be /
    err = ''

    # logger instance
    logger = logging.getLogger(module_name)

    # should there already be a logger instance for the module, return that instance instead
    if (logger.hasHandlers()):
        return logger

    try:
        logger.setLevel(log_level)
    except ValueError:
        logger.setLevel(logging.DEBUG)
        err = f"Unknown logging level: {log_level}, using 'logging.DEBUG instead'"
        log_level = logging.DEBUG

    # formatters to be used by file_handler and console_handler
    file_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(name)s: %(message)s')
    console_formatter = logging.Formatter('%(levelname)s - %(name)s: %(message)s')

    # handlers
    file_handler = logging.FileHandler(log_filepath, mode='a') 
    file_handler.setLevel(log_level)
    file_handler.setFormatter(file_formatter)
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.ERROR) # print to console only ERROR and CRITICAL
    console_handler.setFormatter(console_formatter)

    logger.addHandler(file_handler)
    logger.addHandler(console_handler)
    #logger.propagate = False # Note: this would prohibit log messages to reach the console

    # check if there was an error when setting the log level
    if len(err) > 0:
        logger.error(err)

    return logger