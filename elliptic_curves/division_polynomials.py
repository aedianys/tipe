#!/usr/bin/python3
# -*- coding: utf-8 -*-


class DivisionPolynomialsList:
    def __init__(self, curve_polynomials):
        self.__curve_polynomials = curve_polynomials
        self.__psi = None

    def __getitem__(self, index):
        index = int(index)
        if index < 0:
            raise IndexError

        self.__generate(index)
        return self.__psi[index]

    def psi(self):
        return self.__psi

    def curve_polynomials(self):
        return self.__curve_polynomials

    def __generate(self, l):
        """
        Generate division polynomials `l` if it
        doesn't already exists in self.__psi, using memoization.

        See https://en.wikipedia.org/wiki/Division_polynomials
        for precision on the formulas.
        """
        if not self.__psi:
            self.__psi = 5 * [None]

            # R = F[x,y] / (y**2 - x**3 - A*x - B)
            R = self.__curve_polynomials
            A, B = R.curve.parameters()
            psi = self.__psi

            psi[0] = R(0, 0)
            psi[1] = R(1, 0)
            psi[2] = R(0, 2)
            psi[3] = R((-(A ** 2), 12 * B, 6 * A, 0, 3), 0)
            psi[4] = R(
                0,
                (
                    -4 * (8 * (B ** 2) + A ** 3),
                    -16 * A * B,
                    -20 * (A ** 2),
                    80 * B,
                    20 * A,
                    0,
                    4,
                ),
            )

        psi = self.__psi
        y = self.__curve_polynomials.y()

        # psi[l] is :
        # - a polynomial only in x multiplied by y if l is even
        # - a polynomial only in x if l is odd
        if len(self.__psi) < l + 1:
            self.__psi = self.__psi + [None] * (l + 1 - len(self.__psi))

        if self.__psi[l] is None:
            m, r = divmod(l, 2)
            if r:  # l is odd
                self.__psi[l] = (
                    self[m + 2] * self[m] ** 3 - self[m + 1] ** 3 * self[m - 1]
                )
            else:
                if m % 2 == 0:  # m is even
                    self.__psi[l] = (
                        self[m + 2] * self[m - 1] ** 2 - self[m - 2] * self[m + 1] ** 2
                    ) * (self[m].y_factor() // 2)
                else:  # m is odd
                    self.__psi[l] = (
                        y
                        * (self[m] // 2)
                        * (
                            self[m + 2] * self[m - 1].y_factor() ** 2
                            - self[m - 2] * self[m + 1].y_factor() ** 2
                        )
                    )
