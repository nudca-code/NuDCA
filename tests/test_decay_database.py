#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

from nudca.decay_database import DecayDatabase, DecayDatabaseManager, load_decay_database
from nudca.nuclide import NuclideStrError


class TestDecayDatabase(unittest.TestCase):
    """Test cases for the DecayDatabase class."""

    @classmethod
    def setUpClass(cls):
        """Set up test data that will be used by all test methods."""
        # Create test data
        cls.test_nuclides = np.array(['U-238', 'Th-234', 'Pa-234m', 'U-234'])
        cls.test_half_life_data = np.array([
            (4.468e9 * 365.25 * 24 * 3600, 4.468e9, 'yr', '4.468e9 yr'),  # U-238
            (24.1 * 24 * 3600, 24.1, 'd', '24.1 d'),                      # Th-234
            (1.17 * 60, 1.17, 'min', '1.17 min'),                         # Pa-234m
            (2.455e5 * 365.25 * 24 * 3600, 2.455e5, 'yr', '2.455e5 yr')   # U-234
        ])
        cls.test_decay_modes = np.array([
            ['α'], ['β-'], ['β-'], ['α']
        ])
        cls.test_progeny = np.array([
            ['Th-234'], ['Pa-234m'], ['U-234'], ['Th-230']
        ])
        cls.test_branching_ratios = np.array([
            [1.0], [1.0], [1.0], [1.0]
        ])
        cls.test_decay_energies = np.array([
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.27, 0.0, 0.0, 4.27],  # U-238
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.27, 0.0, 0.0, 0.0, 0.0, 0.27],  # Th-234
            [0.0, 0.0, 0.0, 0.0, 0.0, 2.19, 0.0, 0.0, 0.0, 0.0, 2.19],  # Pa-234m
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.86, 0.0, 0.0, 4.86]   # U-234
        ])

        # Create test database
        cls.db = DecayDatabase(
            data_source='test_data',
            nuclides=cls.test_nuclides,
            half_life_data=cls.test_half_life_data,
            decay_modes_data=cls.test_decay_modes,
            progeny_data=cls.test_progeny,
            branching_ratios_data=cls.test_branching_ratios,
            decay_energies_data=cls.test_decay_energies
        )

    def test_init(self):
        """Test database initialization."""
        self.assertEqual(self.db.data_source, 'test_data')
        self.assertEqual(len(self.db.nuclides), 4)
        self.assertEqual(len(self.db.nuclide_index_map), 4)

    def test_validate_nuclide(self):
        """Test nuclide validation."""
        # Test valid nuclide
        self.assertEqual(self.db._validate_nuclide('U-238'), 'U-238')
        
        # Test invalid format
        with self.assertRaises(NuclideStrError):
            self.db._validate_nuclide('invalid')
        
        # Test non-existent nuclide
        with self.assertRaises(NuclideStrError):
            self.db._validate_nuclide('U-235')

    def test_half_life(self):
        """Test half-life retrieval."""
        # Test readable format
        self.assertEqual(self.db.half_life('U-238'), '4.468e9 yr')
        
        # Test seconds format
        expected_seconds = 4.468e9 * 365.25 * 24 * 3600
        self.assertAlmostEqual(self.db.half_life('U-238', 's'), expected_seconds)
        
        # Test default format
        self.assertEqual(self.db.half_life('U-238', 'yr'), 4.468e9)

    def test_progeny(self):
        """Test progeny retrieval."""
        self.assertEqual(self.db.progeny('U-238'), ['Th-234'])
        self.assertEqual(self.db.progeny('Th-234'), ['Pa-234m'])
        self.assertEqual(self.db.progeny('Pa-234m'), ['U-234'])
        self.assertEqual(self.db.progeny('U-234'), ['Th-230'])

    def test_decay_modes(self):
        """Test decay mode retrieval."""
        self.assertEqual(self.db.decay_modes('U-238'), ['α'])
        self.assertEqual(self.db.decay_modes('Th-234'), ['β-'])
        self.assertEqual(self.db.decay_modes('Pa-234m'), ['β-'])
        self.assertEqual(self.db.decay_modes('U-234'), ['α'])

    def test_decay_energy(self):
        """Test decay energy retrieval."""
        # Test alpha decay energy
        self.assertAlmostEqual(self.db.decay_energy('U-238', 'Alpha'), 4.27)
        self.assertAlmostEqual(self.db.decay_energy('U-234', 'Alpha'), 4.86)
        
        # Test beta decay energy
        self.assertAlmostEqual(self.db.decay_energy('Th-234', 'Beta_Minus'), 0.27)
        self.assertAlmostEqual(self.db.decay_energy('Pa-234m', 'Beta_Minus'), 2.19)
        
        # Test invalid energy type
        with self.assertRaises(ValueError):
            self.db.decay_energy('U-238', 'Invalid')

    def test_decay_energy_shortcuts(self):
        """Test decay energy shortcut methods."""
        self.assertAlmostEqual(self.db.decay_energy_EM('U-238'), 0.0)
        self.assertAlmostEqual(self.db.decay_energy_LP('U-238'), 0.0)
        self.assertAlmostEqual(self.db.decay_energy_HP('U-238'), 0.0)
        self.assertAlmostEqual(self.db.decay_energy_neutrino('U-238'), 0.0)

    def test_plot_nuclear_chart(self):
        """Test nuclear chart plotting."""
        fig = self.db.plot_nuclear_chart(
            min_lifetime=1.0,
            max_lifetime=1e10,
            figsize=(10, 8)
        )
        self.assertIsInstance(fig, plt.Figure)
        plt.close(fig)

    def test_equality(self):
        """Test database equality comparison."""
        # Create identical database
        db2 = DecayDatabase(
            data_source='test_data',
            nuclides=self.test_nuclides,
            half_life_data=self.test_half_life_data,
            decay_modes_data=self.test_decay_modes,
            progeny_data=self.test_progeny,
            branching_ratios_data=self.test_branching_ratios,
            decay_energies_data=self.test_decay_energies
        )
        
        # Create different database
        db3 = DecayDatabase(
            data_source='test_data',
            nuclides=np.array(['U-238']),
            half_life_data=self.test_half_life_data[:1],
            decay_modes_data=self.test_decay_modes[:1],
            progeny_data=self.test_progeny[:1],
            branching_ratios_data=self.test_branching_ratios[:1],
            decay_energies_data=self.test_decay_energies[:1]
        )
        
        self.assertEqual(self.db, db2)
        self.assertNotEqual(self.db, db3)


class TestDecayDatabaseManager(unittest.TestCase):
    """Test cases for the DecayDatabaseManager class."""

    def test_init(self):
        """Test manager initialization."""
        manager = DecayDatabaseManager('test_data')
        self.assertEqual(manager.data_source, 'test_data')

    def test_constants(self):
        """Test constant values."""
        self.assertEqual(DecayDatabaseManager.RADIONUCLIDE_LABEL, 'Radionuclide')
        self.assertEqual(DecayDatabaseManager.MASS_NUMBER_LABEL, 'A')
        self.assertEqual(DecayDatabaseManager.ATOMIC_NUMBER_LABEL, 'Z')


if __name__ == '__main__':
    unittest.main()



