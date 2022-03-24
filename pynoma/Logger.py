import sys
import logging


class Logger:

    logging.basicConfig(stream=sys.stdout, level=logging.INFO)
    handler = logging.getLogger("pynoma")

    @classmethod
    def no_variants_found(self):
        log = "No variants found."
        Logger.handler.info(log)
        return

    @classmethod
    def no_gene_found_with_given_name(self, gene):
        log = f"No gene found with given name: {gene}."
        Logger.handler.info(log)
        return

    @classmethod
    def no_variants_found_for_given_transcript(self, transcript):
        log = f"No variants found for given transcript: {transcript}."
        Logger.handler.info(log)
        return

    @classmethod
    def variant_not_found(self, variant):
        log = f"Variant not found: {variant}."
        Logger.handler.info(log)
        return

    @classmethod
    def batch_searching(self, i, total):
        log = f"Batch searching... {i}/{total}"
        Logger.handler.info(log)
        return
