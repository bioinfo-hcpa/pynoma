import sys
import logging


class Logger:

    logging.basicConfig(stream=sys.stdout, level=logging.INFO)
    handler = logging.getLogger("pynoma")

    @classmethod
    def no_variants_found(cls):
        log = "No variants found."
        Logger.handler.info(log)
        return

    @classmethod
    def no_gene_found_with_given_name(cls, gene):
        log = f"No gene found with given name: {gene}."
        Logger.handler.info(log)
        return

    @classmethod
    def no_variants_found_for_given_transcript(cls, transcript):
        log = f"No variants found for given transcript: {transcript}."
        Logger.handler.info(log)
        return

    @classmethod
    def variant_not_found(cls, variant):
        log = f"Variant not found: {variant}."
        Logger.handler.info(log)
        return

    @classmethod
    def batch_searching(cls, i, total):
        log = f"Batch searching... {i}/{total}"
        Logger.handler.info(log)
        return

    @classmethod
    def request_failed(cls, response):
        log = f"Request failed: {response}."
        Logger.handler.info(log)
        return
    
    @classmethod
    def too_many_requests_error(cls, retry):
        log = "gnomAD is complaining about too many requests. "
        if retry:
            log += f"Retrying in {retry} seconds..."
        else:
            log += "Please try again later. Aborting..."
        Logger.handler.warning(log)
        return
