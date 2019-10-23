#!/usr/bin/env python3
"""Imputation challenge ranking script

Author:
    Jin Lee (leepc12@gmail.com)    
"""

import json
import copy

RANK_PER_CELL_ASSAY = {
    "C19M22": [
        3379312,
        3393574,
        3393756,
        3393847,
        3393457,
        100,
        3330254,
        3388185,
        3393458,
        3393860,
        3393861,
        3389318,
        3391272,
        3393417,
        3393851,
        3379072,
        3393418,
        3393817,
        3393606,
        3386902,
        0,
        3393579,
        3393580,
        3393128,
        3344979
    ],
    "C19M16": [
        3388185,
        3393458,
        3393574,
        3330254,
        3393457,
        3393847,
        3393756,
        3393606,
        3393817,
        3379072,
        3393418,
        3379312,
        3393417,
        3393580,
        3393860,
        3393861,
        100,
        3393128,
        3393579,
        3389318,
        3386902,
        3391272,
        0,
        3393851,
        3344979
    ],
    "C38M01": [
        3393574,
        3393847,
        3330254,
        3393417,
        100,
        3379072,
        3393756,
        3393418,
        3393457,
        3391272,
        3393817,
        3393860,
        3393861,
        3393458,
        3379312,
        3393580,
        3388185,
        3389318,
        0,
        3386902,
        3393606,
        3393579,
        3393128,
        3393851,
        3344979
    ],
    "C51M16": [
        3330254,
        3393851,
        3393574,
        3393417,
        3389318,
        3393847,
        3388185,
        3393606,
        100,
        3391272,
        3393458,
        3393756,
        3393457,
        3393418,
        3379072,
        3393817,
        3393860,
        3393861,
        3379312,
        3393580,
        3393128,
        3386902,
        0,
        3393579,
        3344979
    ],
    "C06M16": [
        3330254,
        3393574,
        3393847,
        3393756,
        3388185,
        100,
        3393457,
        3393417,
        3393851,
        3391272,
        3393606,
        0,
        3379312,
        3389318,
        3393458,
        3393817,
        3379072,
        3393418,
        3393860,
        3393861,
        3393128,
        3393580,
        3386902,
        3393579,
        3344979
    ],
    "C40M22": [
        100,
        3389318,
        3391272,
        3393417,
        3393851,
        3379072,
        3393847,
        3393574,
        3379312,
        3393458,
        3393860,
        3393861,
        3393457,
        3393817,
        3388185,
        3393606,
        3330254,
        3393756,
        3393579,
        3386902,
        3393418,
        0,
        3393580,
        3393128,
        3344979
    ],
    "C06M17": [
        3393580,
        3393457,
        3393817,
        3393756,
        3393417,
        3393418,
        3393458,
        3393860,
        3393861,
        3388185,
        3393574,
        3393847,
        3379072,
        100,
        3330254,
        3393606,
        3379312,
        3386902,
        3393579,
        3389318,
        3344979,
        0,
        3393851,
        3393128,
        3391272
    ],
    "C31M01": [
        3379072,
        3393860,
        3393861,
        3393458,
        100,
        3330254,
        3393417,
        3393418,
        3393574,
        3393847,
        3393851,
        3379312,
        3393756,
        0,
        3389318,
        3393580,
        3393457,
        3344979,
        3386902,
        3393606,
        3393579,
        3393817,
        3388185,
        3391272,
        3393128
    ],
    "C07M29": [
        3393817,
        3379072,
        3393457,
        3393606,
        3393756,
        3393458,
        3393860,
        3393861,
        100,
        3388185,
        3330254,
        3393417,
        3393580,
        3393418,
        3393574,
        3393847,
        3379312,
        3389318,
        3393579,
        0,
        3386902,
        3344979,
        3393851,
        3393128,
        3391272
    ],
    "C19M17": [
        3393457,
        3393606,
        3379072,
        3393817,
        3393458,
        3393860,
        3393861,
        3393756,
        3389318,
        3393580,
        3393574,
        3386902,
        3393847,
        3344979,
        3379312,
        3388185,
        3330254,
        3393417,
        3393418,
        3393579,
        100,
        0,
        3393128,
        3391272,
        3393851
    ],
    "C38M20": [
        3393817,
        3330254,
        3393417,
        3393606,
        3379312,
        3388185,
        3389318,
        3393580,
        3393457,
        3393847,
        3393418,
        3393574,
        3379072,
        3393756,
        3393458,
        3386902,
        3344979,
        3393579,
        100,
        3393128,
        3393860,
        3393851,
        0,
        3393861,
        3391272
    ],
    "C05M29": [
        3393458,
        3393860,
        3393861,
        3379072,
        3393817,
        100,
        3393606,
        3330254,
        3393457,
        3393417,
        3393580,
        3393756,
        3379312,
        3389318,
        3386902,
        3393579,
        3393418,
        3388185,
        3344979,
        0,
        3393574,
        3393847,
        3393128,
        3393851,
        3391272
    ],
    "C28M17": [
        3393756,
        3393457,
        3393574,
        3393847,
        3389318,
        3393817,
        3393580,
        3388185,
        3393458,
        3393860,
        3393861,
        3393606,
        3379072,
        3386902,
        3379312,
        3344979,
        3330254,
        3393417,
        3393579,
        100,
        3393418,
        0,
        3393128,
        3391272,
        3393851
    ],
    "C22M17": [
        3379072,
        3393817,
        3393458,
        3393860,
        3393861,
        3393606,
        3393580,
        3393756,
        3393457,
        3388185,
        3393847,
        3379312,
        3393574,
        3393418,
        3393417,
        3393579,
        3330254,
        3386902,
        3344979,
        3389318,
        100,
        0,
        3393128,
        3391272,
        3393851
    ],
    "C12M01": [
        3393457,
        100,
        3391272,
        3393847,
        3393574,
        3330254,
        3393417,
        3379072,
        3393756,
        3393860,
        3393861,
        3393458,
        3388185,
        3393851,
        3393418,
        3389318,
        3379312,
        3386902,
        3344979,
        3393580,
        0,
        3393817,
        3393579,
        3393128,
        3393606
    ],
    "C38M22": [
        3393457,
        3393756,
        3391272,
        3393847,
        100,
        3393417,
        3379312,
        3330254,
        3393574,
        3388185,
        0,
        3393418,
        3393458,
        3393860,
        3393861,
        3393606,
        3379072,
        3393817,
        3393851,
        3386902,
        3389318,
        3393579,
        3393580,
        3393128,
        3344979
    ],
    "C19M18": [
        3393756,
        3393580,
        3393457,
        100,
        3393417,
        3330254,
        3388185,
        3393860,
        3393861,
        3393458,
        3393847,
        3393574,
        3393817,
        3393418,
        3379072,
        3393606,
        3389318,
        3379312,
        3393579,
        3386902,
        0,
        3393128,
        3344979,
        3391272,
        3393851
    ],
    "C31M16": [
        3330254,
        3393417,
        100,
        3393574,
        3388185,
        3393847,
        3379312,
        3393860,
        3393861,
        3393851,
        3393458,
        3389318,
        3391272,
        0,
        3379072,
        3393580,
        3393817,
        3393418,
        3393457,
        3393579,
        3393606,
        3393128,
        3386902,
        3393756,
        3344979
    ],
    "C28M18": [
        100,
        3388185,
        3393580,
        3393417,
        3330254,
        3393860,
        3393861,
        3393458,
        3393847,
        3393756,
        3393817,
        3393574,
        3379312,
        3393418,
        3379072,
        3393457,
        3393606,
        3389318,
        3386902,
        3344979,
        0,
        3393128,
        3393579,
        3391272,
        3393851
    ],
    "C28M29": [
        3393817,
        3393458,
        3393860,
        3393861,
        3393457,
        100,
        3393580,
        3393756,
        3393574,
        3379072,
        3330254,
        3393606,
        3389318,
        3393847,
        3393417,
        3393418,
        3388185,
        3386902,
        0,
        3393579,
        3393128,
        3379312,
        3393851,
        3344979,
        3391272
    ],
    "C28M22": [
        3379312,
        3393574,
        3393847,
        3330254,
        100,
        3393457,
        3393756,
        3389318,
        3388185,
        3379072,
        3393606,
        3393417,
        3393458,
        3393860,
        3393861,
        3391272,
        3393817,
        3386902,
        0,
        3393418,
        3393851,
        3393579,
        3393580,
        3393128,
        3344979
    ],
    "C38M29": [
        3393457,
        3393756,
        3389318,
        3393817,
        3393458,
        3393860,
        3393861,
        3330254,
        3393417,
        3388185,
        3393418,
        3393580,
        100,
        3393574,
        3379072,
        3393606,
        3393847,
        3386902,
        3393579,
        3379312,
        3344979,
        3393128,
        3391272,
        0,
        3393851
    ],
    "C51M29": [
        3393817,
        3393458,
        3393860,
        3393861,
        100,
        3330254,
        3393417,
        3379072,
        3393606,
        3393580,
        3393418,
        3388185,
        3393457,
        3393756,
        3393128,
        3379312,
        3393574,
        3386902,
        3393847,
        3393579,
        3389318,
        3344979,
        0,
        3391272,
        3393851
    ],
    "C51M20": [
        3393580,
        3393417,
        3393418,
        3389318,
        3388185,
        3393574,
        3393847,
        3330254,
        3393817,
        3393756,
        3379312,
        3393457,
        3393606,
        3393458,
        3379072,
        3393860,
        0,
        3393861,
        3386902,
        3344979,
        3393128,
        100,
        3393579,
        3393851,
        3391272
    ],
    "C40M16": [
        3393417,
        100,
        3388185,
        3393756,
        3393457,
        3393574,
        3393851,
        3391272,
        3393847,
        3389318,
        3379072,
        3393606,
        3393418,
        3393860,
        3393861,
        3393458,
        3393580,
        3393817,
        3393128,
        0,
        3379312,
        3386902,
        3393579,
        3330254,
        3344979
    ],
    "C39M29": [
        3393457,
        3393756,
        3393847,
        3388185,
        3393574,
        3393817,
        3379072,
        3393606,
        3389318,
        100,
        3393458,
        3393860,
        3393861,
        3344979,
        3330254,
        3393580,
        3393417,
        3393418,
        3379312,
        3386902,
        3393851,
        0,
        3393579,
        3393128,
        3391272
    ],
    "C38M02": [
        100,
        3393574,
        3393817,
        3393847,
        3379072,
        3393458,
        3393860,
        3393417,
        3393861,
        0,
        3330254,
        3389318,
        3393606,
        3379312,
        3393418,
        3393756,
        3388185,
        3393580,
        3386902,
        3393851,
        3393457,
        3393579,
        3393128,
        3344979,
        3391272
    ],
    "C22M16": [
        3330254,
        3391272,
        3393574,
        3393847,
        3393851,
        3388185,
        100,
        3393417,
        3393756,
        3393457,
        3393458,
        3389318,
        3379312,
        3379072,
        3393606,
        3393817,
        3393418,
        0,
        3393579,
        3393860,
        3393861,
        3393580,
        3386902,
        3393128,
        3344979
    ],
    "C40M18": [
        100,
        3393417,
        3393580,
        3393457,
        3393574,
        3393847,
        3393860,
        3393861,
        3393458,
        3330254,
        3388185,
        3393128,
        3393418,
        3393756,
        3393817,
        3379072,
        3393606,
        3379312,
        3389318,
        0,
        3386902,
        3393579,
        3344979,
        3393851,
        3391272
    ],
    "C39M20": [
        3389318,
        3330254,
        3388185,
        3393847,
        3393574,
        3393580,
        3393817,
        3393756,
        3393417,
        3393606,
        3393457,
        3379312,
        3379072,
        3393418,
        3393458,
        3386902,
        3344979,
        3393579,
        3393860,
        3393128,
        0,
        3393861,
        100,
        3393851,
        3391272
    ],
    "C19M20": [
        3389318,
        3393817,
        3393580,
        3393606,
        3379312,
        3393574,
        3393847,
        3379072,
        3388185,
        3393756,
        3393418,
        3344979,
        3393457,
        3393458,
        3330254,
        3393417,
        3386902,
        3393579,
        0,
        3393861,
        3393860,
        3393128,
        100,
        3393851,
        3391272
    ],
    "C07M20": [
        3393574,
        3393847,
        3388185,
        3330254,
        3393817,
        3393580,
        3393418,
        3379072,
        3393458,
        3393417,
        3393606,
        100,
        3379312,
        3393457,
        3393579,
        3389318,
        0,
        3393861,
        3393756,
        3393860,
        3393128,
        3386902,
        3344979,
        3393851,
        3391272
    ],
    "C19M29": [
        3393817,
        3393458,
        3393860,
        3393861,
        3393457,
        3379312,
        3393756,
        3330254,
        100,
        3393417,
        3379072,
        3393606,
        3389318,
        3393418,
        3388185,
        3393580,
        3386902,
        3393574,
        3344979,
        3393847,
        3393851,
        3393579,
        3393128,
        0,
        3391272
    ],
    "C51M18": [
        100,
        3393457,
        3393417,
        3393756,
        3330254,
        3393580,
        3393418,
        3379072,
        3393606,
        3393817,
        3393860,
        3393861,
        3393458,
        3393128,
        3388185,
        3393574,
        3389318,
        3393847,
        3393579,
        3386902,
        3379312,
        3344979,
        0,
        3393851,
        3391272
    ],
    "C38M17": [
        3330254,
        3393417,
        100,
        3393418,
        3379072,
        3393606,
        3393458,
        3393860,
        3393861,
        3388185,
        3393756,
        3393457,
        3393817,
        3393128,
        3393580,
        3393579,
        3393574,
        3393847,
        3344979,
        3386902,
        3379312,
        0,
        3389318,
        3393851,
        3391272
    ],
    "C39M22": [
        100,
        3389318,
        3393417,
        3379312,
        3393851,
        3391272,
        3393458,
        3393860,
        3393861,
        3393606,
        3330254,
        3379072,
        3393574,
        3393847,
        3393418,
        3388185,
        3393817,
        3393457,
        0,
        3393756,
        3386902,
        3393128,
        3393580,
        3393579,
        3344979
    ],
    "C31M02": [
        100,
        3393458,
        3393851,
        3379072,
        3393860,
        3393574,
        3393861,
        3393847,
        3393817,
        3393417,
        3330254,
        3379312,
        3389318,
        3393606,
        0,
        3393418,
        3388185,
        3393756,
        3344979,
        3386902,
        3393457,
        3391272,
        3393580,
        3393579,
        3393128
    ],
    "C05M18": [
        3393817,
        3393606,
        3393860,
        3393861,
        3393458,
        3393847,
        3393457,
        3379072,
        3393574,
        100,
        3393756,
        3330254,
        3393580,
        3388185,
        3393417,
        3393418,
        3379312,
        3393579,
        3389318,
        3386902,
        0,
        3393128,
        3344979,
        3393851,
        3391272
    ],
    "C39M18": [
        3393457,
        3393847,
        3388185,
        3393574,
        3393756,
        3393817,
        3379312,
        3393606,
        3393860,
        3393861,
        3393458,
        3379072,
        100,
        3393417,
        3393580,
        3393418,
        3330254,
        3389318,
        3393851,
        3386902,
        3393579,
        3344979,
        0,
        3393128,
        3391272
    ],
    "C51M17": [
        3379072,
        3393418,
        3393458,
        3393860,
        3393861,
        3393417,
        3393457,
        3393817,
        3330254,
        3393580,
        3393606,
        3393756,
        100,
        3388185,
        3393574,
        3379312,
        3393847,
        3386902,
        3393579,
        3344979,
        3389318,
        3393128,
        3393851,
        0,
        3391272
    ],
    "C12M02": [
        3389318,
        3393847,
        3393574,
        3330254,
        3393417,
        3388185,
        3379072,
        3393458,
        3393756,
        3393817,
        3391272,
        3393861,
        3393457,
        100,
        3393606,
        3393860,
        3379312,
        3393851,
        3393418,
        3386902,
        3344979,
        3393580,
        0,
        3393579,
        3393128
    ],
    "C38M18": [
        3393574,
        3393847,
        3388185,
        3393580,
        100,
        3393417,
        3393756,
        3330254,
        3393817,
        0,
        3389318,
        3393579,
        3393128,
        3393418,
        3393457,
        3393606,
        3379072,
        3393860,
        3393861,
        3393458,
        3379312,
        3386902,
        3344979,
        3393851
    ],
    "C40M29": [
        3393817,
        3393458,
        3393860,
        3393861,
        3379072,
        3393756,
        3393606,
        3388185,
        3393457,
        3393580,
        3389318,
        100,
        3330254,
        3393417,
        3393847,
        3386902,
        3393418,
        3393574,
        3393579,
        3379312,
        3344979,
        0,
        3393851,
        3391272,
        3393128
    ],
    "C06M18": [
        3393580,
        3393817,
        100,
        3393417,
        3330254,
        3379072,
        3393606,
        3393860,
        3393861,
        3393458,
        3393418,
        3379312,
        3393457,
        3388185,
        3389318,
        3393756,
        3393579,
        3393574,
        3386902,
        3393847,
        3344979,
        3393128,
        0,
        3393851,
        3391272
    ],
    "C05M17": [
        3393417,
        100,
        3330254,
        3393418,
        3393580,
        3379312,
        3393817,
        3393756,
        3393574,
        3393847,
        3393458,
        3393860,
        3393861,
        3393457,
        3388185,
        3379072,
        3393128,
        0,
        3393606,
        3391272,
        3344979,
        3386902,
        3393579,
        3389318,
        3393851
    ],
    "C31M29": [
        3393606,
        3393817,
        3393580,
        3388185,
        3393458,
        3393860,
        3393861,
        3379072,
        3393574,
        3393457,
        3389318,
        3393418,
        3393847,
        3393417,
        3330254,
        100,
        3393756,
        3379312,
        3386902,
        3393579,
        3393128,
        3344979,
        0,
        3393851,
        3391272
    ],
    "C05M20": [
        3389318,
        3393580,
        3393817,
        3393458,
        3393417,
        3393418,
        3344979,
        3379072,
        3388185,
        3330254,
        3393606,
        3386902,
        3379312,
        3393847,
        3393457,
        3393756,
        3393574,
        3393860,
        3393128,
        0,
        3393861,
        3393579,
        100,
        3393851,
        3391272
    ],
    "C40M20": [
        3393417,
        3389318,
        3393580,
        3388185,
        3393574,
        3393847,
        3393128,
        3393457,
        3393756,
        100,
        3393860,
        3330254,
        3379312,
        3393851,
        3393418,
        3393579,
        0,
        3393861,
        3393458,
        3393817,
        3393606,
        3379072,
        3391272,
        3386902,
        3344979
    ],
    "C40M17": [
        3388185,
        3393574,
        3393756,
        3393580,
        3393847,
        3393457,
        3393817,
        3393417,
        3379312,
        3393458,
        3393860,
        3393861,
        3393606,
        100,
        3330254,
        3379072,
        0,
        3389318,
        3393851,
        3393418,
        3386902,
        3393579,
        3344979,
        3391272,
        3393128
    ],
    "C39M17": [
        3393756,
        3393457,
        3389318,
        3386902,
        3344979,
        3388185,
        3393847,
        3393574,
        3393817,
        3393579,
        3379312,
        3393606,
        3379072,
        3393458,
        3393860,
        3393861,
        3393580,
        100,
        3393417,
        0,
        3330254,
        3393418,
        3393851,
        3393128,
        3391272
    ],
    "C39M16": [
        3393417,
        3393851,
        100,
        3391272,
        3393847,
        3393574,
        3389318,
        3388185,
        3393128,
        3393418,
        3393457,
        3393606,
        3379072,
        3393756,
        3393458,
        3393817,
        3393580,
        3379312,
        3386902,
        3393860,
        3393861,
        3330254,
        0,
        3393579,
        3344979
    ]
}

TEAM_NAME = {
    -1:      "Blind",
    0:       "Avocado",
    100:     "Average",
    3330254: "Hongyang Li and Yuanfang Guan",
    3393128: "Aug2019Imputation",
    3388185: "LiPingChun",
    3393417: "Hongyang Li and Yuanfang Guan v1",
    3393418: "Hongyang Li and Yuanfang Guan v2",
    3379072: "BrokenNodes",
    3391272: "UIOWA Michaelson Lab",
    3379312: "CostaLab",
    3389318: "Song Lab",
    3393574: "Lavawizard",
    3393457: "imp",
    3393756: "imp1",
    3393606: "BrokenNodes_v2",
    3393817: "BrokenNodes_v3",
    3344979: "NittanyLions",
    3393458: "CUWA",
    3393579: "Song Lab 2",
    3393580: "Song Lab 3",
    3393847: "Guacamole",
    3393860: "CUImpute1",
    3393861: "ICU",
    3393851: "NittanyLions2",
    3386902: "KKT-ENCODE-Impute",
}

TEAM_COLOR = {
    -1:      (0, 0, 0),
    0:       (128, 0, 0),
    100:     (170, 110, 40),
    3330254: (128, 128, 0),
    3393128: (0, 128, 128),
    3388185: (0, 0, 128),
    3393417: (255, 255, 0),
    3393418: (230, 25, 75),
    3379072: (245, 130, 48),
    3391272: (255, 225, 25),
    3379312: (210, 245, 60),
    3389318: (60, 180, 75),
    3393574: (70, 240, 240),
    3393457: (0, 130, 200),
    3393756: (145, 30, 180),
    3393606: (240, 50, 230),
    3393817: (128, 128, 128),
    3344979: (250, 190, 190),
    3393458: (255, 215, 180),
    3393579: (255, 250, 200),
    3393580: (170, 255, 195),
    3393847: (230, 190, 255),
    3393860: (0, 0, 255),
    3393861: (0, 255, 0),
    3393851: (255, 0, 0),
    3386902: (0, 255, 128),
}

URL_ROOT = {
    -1:     "http://mitra.stanford.edu/kundaje/ic/blind",
    0:      "http://mitra.stanford.edu/kundaje/ic/avocado",
    100:    "http://mitra.stanford.edu/kundaje/ic/average",
}

for k, v in TEAM_NAME.items():
    if k < 1000:
        continue
    URL_ROOT[k] = "http://mitra.stanford.edu/kundaje/ic/round2/" + str(k)

#{"type":"coordinate_override","coord":"chr9,36329955,chr9,37537411"},
for cell_assay, team_ids in RANK_PER_CELL_ASSAY.items():
    hub_url = "http://mitra.stanford.edu/kundaje/ic/track_hubs/"+cell_assay+".txt"
    print("http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&hgct_customText=" + hub_url)
    with open(cell_assay + '.txt', 'w') as fp:
        STR='''
hub ENCODE IC BWs for {cell_assay}
shortLabel ic_{cell_assay}
longLabel ENCODE Imputation Challenge signal tracks (bigwigs) for {cell_assay}
'''.format(cell_assay=cell_assay)
        STR=''
        for i, k in enumerate([-1] + team_ids):
            if cell_assay=='C38M18' and k==3391272:
                # missing submission of team 3391272
                continue
            r, g, b = TEAM_COLOR[k]
            STR+='track type=bigWig name="{name}" priority={priority} smoothingWindow=3 color={r},{g},{b} autoScale=off viewLimits=0:40 visibility=full windowingFunction=maximum bigDataUrl={url}\n'.format(
                name='{name} ({id}, {c})'.format(id=k, name=TEAM_NAME[k], c=cell_assay),
                url='{root}/{cell_assay}.bigwig'.format(root=URL_ROOT[k], cell_assay=cell_assay),
                r=r,
                g=g,
                b=b,
                priority=i+1)
        fp.write(STR)
