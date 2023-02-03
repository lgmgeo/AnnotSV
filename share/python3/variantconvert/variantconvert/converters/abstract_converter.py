# -*- coding: utf-8 -*-

import json

from abc import ABC, abstractmethod


class AbstractConverter(ABC):
    def __init__(self, config_filepath):
        self.config_filepath = config_filepath
        with open(config_filepath, "r") as f:
            self.config = json.load(f)

    @abstractmethod
    def convert(self, file, output_path):
        pass
