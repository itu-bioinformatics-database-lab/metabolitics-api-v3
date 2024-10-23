import os
import json
import gzip

import cobra as cb
import pandas as pd
import escher
from cobra.core import Model, DictList, Reaction, Metabolite


class DataReader(object):
    def __init__(self):
        self.path = '../dataset/disease'
        self.y_label = 'stage'

    def read_data(self, disease_name, by_stage=False):
        df = pd.read_csv('%s/%s.csv' % (self.path, disease_name), header=0)
        X = df.ix[:, df.columns != self.y_label].to_dict('records')
        y = df[self.y_label].values
        if not by_stage:
            y = ['h' if i == 'h' else disease_name.lower() for i in y]
        return (X, y)

    def read_disease_sample(self, disease_name):
        X, y = self.read_data(disease_name)
        X_new, y_new = [], []
        for label in ['h', disease_name.lower()]:
            x, l = list(filter(lambda x: x[1] == label, zip(X, y)))[0]
            X_new.append(x)
            y_new.append(l)
        return X_new, y_new

    def read_healthy(self, disease_name):
        return list(
            zip(*filter(lambda y: y[1] == 'h',
                        zip(*DataReader().read_data(disease_name)))))

    def read_columns(self, disease_name):
        return pd.read_csv(
            '%s/%s.csv' % (self.path, disease_name), header=0).columns

    def read_all(self):
        disease_names = ['BC', 'HCC']
        (X, y) = (list(), list())
        for i in disease_names:
            (tX, ty) = self.read_data(i)
            X += tX
            y += ty
        return (X, y)

    def read_avg_data(self):
        raise NotImplementedError()

    def read_small_data(self):
        self.path = '../dataset/small-disease'
        return self.read_all()

    def read_escher_map(self, map_name):
        return escher.Builder(
            map_json='../dataset/visualizations/%s_map.json' % map_name)

    def read_categorical_solutions(self):
        raise NotImplemented()

    def read_network_model(self, name='recon3D'):
        return cb.io.load_json_model('../dataset/network/%s.json' % name)

    def read_subsystem_categories(self, name='recon'):
        path = '../dataset/subsystem-categories/%s.json' % name
        with open(path) as f:
            return {k: set(v) for k, v in json.load(f).items()}

    def create_example_model(self):
        model = Model('example_model')

        rs = (r1, r2, r3) = (Reaction('R1'), Reaction('R2'), Reaction('R3'))
        for r in rs:
            r.lower_bound = 0.
            r.upper_bound = 1000.
            r.objective_coefficient = 0.

        ACP_c = Metabolite(
            'ACP_c',
            formula='C11H21N2O7PRS',
            name='acyl-carrier-protein',
            compartment='c')

        r1.add_metabolites({ACP_c: 1.0})
        r2.add_metabolites({ACP_c: 1.0})
        r3.add_metabolites({ACP_c: -1.0})

        model.add_reactions(rs)

        return model

    def create_example2_model(self):
        model = Model('example2_model')

        rs = (r1, r2, r3, r4,
              r5) = (Reaction('R1'), Reaction('R2'), Reaction('R3'),
                     Reaction('R4'), Reaction('R5'))
        for r in rs:
            r.lower_bound = 0.
            r.upper_bound = 1000.
            r.objective_coefficient = 0.

        ms = (m1, m2, m3) = (Metabolite('M1'), Metabolite('M2'),
                             Metabolite('M3'))

        r1.add_metabolites({m1: 1.0})
        r2.add_metabolites({m1: -1.0, m2: 1.0})
        r3.add_metabolites({m2: -1.0})
        r4.add_metabolites({m1: -1.0, m3: 1.0})
        r5.add_metabolites({m3: -1.0})

        model.add_reactions(rs)

        return model

    def read_hmdb_diseases(self):
        return self.read_json('%s/hmdb_disease_measurements.json' % self.path)

    def read_json(self, path, gz=False):
        with gzip.open('%s.gz' % path, 'rt') if gz else open(path) as f:
            return json.load(f)

    def read_solution(self, filename):
        solution = self.read_json('../dataset/solutions/%s.json' % filename)
        return [solution.values(), solution.keys()]

    def read_analyze_solution(self, filename, gz=True):
        path = '../dataset/solutions/%s.json' % filename
        return list(zip(*self.read_json(path, gz=gz)))
