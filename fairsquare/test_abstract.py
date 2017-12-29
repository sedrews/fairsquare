import unittest
from abstract import *
from fractions import Fraction

class TestAbstract(unittest.TestCase):
    
    # Test Tristate operations

    def test_tristate_algebra_1(self):
        t,f,u = Tristate.TRUE,Tristate.FALSE,Tristate.UNKNOWN

        self.assertEqual(t & t, t)
        self.assertEqual(t & f, f)
        self.assertEqual(t & u, u)
        self.assertEqual(f & t, f)
        self.assertEqual(f & f, f)
        self.assertEqual(f & u, f)
        self.assertEqual(u & t, u)
        self.assertEqual(u & u, u)
        self.assertEqual(u & f, f)

        self.assertEqual(t | t, t)
        self.assertEqual(t | f, t)
        self.assertEqual(t | u, t)
        self.assertEqual(f | t, t)
        self.assertEqual(f | f, f)
        self.assertEqual(f | u, u)
        self.assertEqual(u | t, t)
        self.assertEqual(u | f, u)
        self.assertEqual(u | u, u)

        self.assertEqual(-t, f)
        self.assertEqual(-f, t)
        self.assertEqual(-u, u)

    def test_tristate_algebra_2(self):
        t,f,u = Tristate.TRUE,Tristate.FALSE,Tristate.UNKNOWN
        
        # Test __rand__
        self.assertEqual(True & t, t)
        self.assertEqual(True & f, f)
        self.assertEqual(True & u, u)
        self.assertEqual(False & t, f)
        self.assertEqual(False & f, f)
        self.assertEqual(False & u, f)

        # Test __ror__
        self.assertEqual(True | t, t)
        self.assertEqual(True | f, t)
        self.assertEqual(True | u, t)
        self.assertEqual(False | t, t)
        self.assertEqual(False | f, f)
        self.assertEqual(False | u, u)

    def test_tristate_to_bool(self):
        t,f,u = Tristate.TRUE,Tristate.FALSE,Tristate.UNKNOWN
        self.assertEqual(bool(t), True)
        self.assertEqual(bool(f), False)
        self.assertRaises(ValueError, Tristate.__bool__, u)

    def test_tristate_algebra_typeerror(self):
        self.assertRaises(TypeError, lambda x,y : x & y, Tristate.TRUE, 5)
        self.assertRaises(TypeError, lambda x,y : x & y, 'foo', Tristate.FALSE)
        self.assertRaises(TypeError, lambda x,y : x | y, Tristate.UNKNOWN, 1.1)
        self.assertRaises(TypeError, lambda x,y : x | y, ['bar', 'baz'], Tristate.TRUE)
        self.assertRaises(AssertionError, Tristate.__neg__, lambda x : x)
        
    # Test Interval auxiliary methods

    def test_interval_normal_form(self):
        self.assertEqual(Interval._normalize([(1,2), (0,4)]), [(0,4)])
        self.assertEqual(Interval._normalize([(0,2), (1,1), (2,3), (4.5,4.5)]),
                         [(0,3), (4.5,4.5)])

    def test_interval_init(self):
        self.assertTrue(len(Interval((0,4),(1,2))._pairs) == 1)
        self.assertRaises(AssertionError, Interval)
        self.assertRaises(AssertionError, Interval, (0,1), (3,2))

    def test_interval_aux(self):
        # Test __contains__
        self.assertTrue(0 in Interval((0,0)))
        self.assertTrue(1 in Interval((-1,0), (0.5, 3)))
        self.assertFalse(0 in Interval((5,10)))
        self.assertRaises(ValueError, lambda x : x in Interval((7,8)), "hello")

        # Test pairs property
        i = Interval((0,0))
        i.pairs = [(0,2), (1,3)]
        self.assertEqual(i.pairs, [(0,3)])

    # Test Interval comparisons

    def test_interval_concrete_comparisons(self):
        a = Interval((0,1), (3,4))
        b = Interval((4,4))
        c = Interval((6,6), (7,8))
        d = Interval((3,6.5))

        t,f,u = Tristate.TRUE,Tristate.FALSE,Tristate.UNKNOWN

        self.assertEqual(a < a, u)
        self.assertEqual(a < b, u)
        self.assertEqual(a < c, t)
        self.assertEqual(a < d, u)
        self.assertEqual(b < a, f)
        self.assertEqual(b < b, f)
        self.assertEqual(b < c, t)
        self.assertEqual(b < d, u)
        self.assertEqual(c < a, f)
        self.assertEqual(c < b, f)
        self.assertEqual(c < c, u)
        self.assertEqual(c < d, u)
        self.assertEqual(d < a, u)
        self.assertEqual(d < b, u)
        self.assertEqual(d < c, u)
        self.assertEqual(d < d, u)

        self.assertEqual(a <= a, u)
        self.assertEqual(a <= b, t)
        self.assertEqual(a <= c, t)
        self.assertEqual(a <= d, u)
        self.assertEqual(b <= a, u)
        self.assertEqual(b <= b, t)
        self.assertEqual(b <= c, t)
        self.assertEqual(b <= d, u)
        self.assertEqual(c <= a, f)
        self.assertEqual(c <= b, f)
        self.assertEqual(c <= c, u)
        self.assertEqual(c <= d, u)
        self.assertEqual(d <= a, u)
        self.assertEqual(d <= b, u)
        self.assertEqual(d <= c, u)
        self.assertEqual(d <= d, u)

        self.assertEqual(a == a, u)
        self.assertEqual(a == b, u)
        self.assertEqual(a == c, f)
        self.assertEqual(a == d, u)
        self.assertEqual(b == b, t)
        self.assertEqual(b == c, f)
        self.assertEqual(b == d, u)
        self.assertEqual(c == c, u)
        self.assertEqual(c == d, u)
        self.assertEqual(d == d, u)

        self.assertEqual(a >= a, a <= a)
        self.assertEqual(a > a, a < a)
        self.assertEqual(a > b, b < a)
        self.assertEqual(a > c, c < a)
        self.assertEqual(a > d, d < a)
        self.assertEqual(b >= a, a <= b)
        self.assertEqual(b >= b, b <= b)
        self.assertEqual(b > b, b < b)
        self.assertEqual(b > c, c < b)
        self.assertEqual(b > d, d < b)
        self.assertEqual(c >= a, a <= c)
        self.assertEqual(c >= b, b <= c)
        self.assertEqual(c >= c, c <= c)
        self.assertEqual(c > c, c < c)
        self.assertEqual(c > d, d < c)
        self.assertEqual(d >= a, a <= d)
        self.assertEqual(d >= b, b <= d)
        self.assertEqual(d >= c, c <= d)
        self.assertEqual(d >= d, d <= d)
        self.assertEqual(d > d, d < d)

    def test_interval_comparison_typeerror(self):
        self.assertEqual(Interval((0,1)) < 2, Tristate.TRUE)
        self.assertEqual(Interval((0,0)) == 0, Tristate.TRUE)
        self.assertRaises(TypeError, lambda x,y : x >= y, Interval((0,1)), "foo")

    def test_interval_abstract_equality(self):
        i1 = Interval((0,1), (3,5), (4,5), (2,2), (3,6))
        i2 = Interval((0,1), (2,2), (3,6))
        self.assertEqual(i1 == i2, Tristate.UNKNOWN)
        self.assertTrue(i1.pairs == i2.pairs)

    # Test Interval arithmetic

    def test_interval_negation(self):
        i1 = -Interval((0,1))
        i2 = Interval((-1,0))
        self.assertTrue(i1.pairs == i2.pairs)
        i1 = -Interval((-2,3), (4,5))
        i2 = Interval((-5,-4), (-3,2))
        self.assertTrue(i1.pairs == i2.pairs)

    def test_interval_addition(self):
        a = Interval((0,2), (4,5))
        b = Interval((-2, 2))
        c = Interval((-2, 7))
        self.assertTrue((a+b).pairs == c.pairs)
        d = Interval((1,3), (5,6))
        self.assertTrue((a+1).pairs == d.pairs)

    def test_interval_subtraction(self):
        a = Interval((0,2), (4,5))
        b = Interval((1,2))
        c = Interval((-2,1), (2,4))
        self.assertTrue((a-b).pairs == c.pairs)
        d = Interval((-1,1), (3,4))
        self.assertTrue((a-1).pairs == d.pairs)

    def test_interval_multiplication(self):
        a = Interval((1,2))
        b = Interval((-3,4))
        c = Interval((-6,-5))
        self.assertEqual((a*a).pairs, [(1,4)])
        self.assertEqual((a*b).pairs, [(-6,8)])
        self.assertEqual((a*c).pairs, [(-12,-5)])
        self.assertEqual((b*b).pairs, [(-12,16)])
        self.assertEqual((b*c).pairs, [(-24,18)])
        self.assertEqual((c*c).pairs, [(25,36)])

    def test_interval_division(self):
        a = Interval((Fraction(1),Fraction(2)))
        b = Interval((Fraction(-3),Fraction(4)))
        c = Interval((Fraction(-6),Fraction(-5)))
        self.assertEqual((a/a).pairs, [(Fraction(1,2),Fraction(2))])
        self.assertRaises(ZeroDivisionError, lambda x,y : x / y, a, b)
        self.assertEqual((a/c).pairs, [(-Fraction(2,5),-Fraction(1,6))])
        self.assertEqual((b/a).pairs, [(Fraction(-3),Fraction(4))])
        self.assertRaises(ZeroDivisionError, lambda x,y : x / y, b, b)
        self.assertEqual((b/c).pairs, [(-Fraction(4,5),Fraction(3,5))])
        self.assertEqual((c/a).pairs, [(-Fraction(6),-Fraction(5,2))])
        self.assertRaises(ZeroDivisionError, lambda x,y : x / y, c, b)
        self.assertEqual((c/c).pairs, [(Fraction(5,6),Fraction(6,5))])

    def test_interval_arithmetic_r_methods(self):
        # Test the radd rsub etc on numeric types
        i = Interval((1,2))
        self.assertEqual((4+i).pairs, [(5,6)])
        self.assertEqual((4-i).pairs, [(2,3)])
        self.assertEqual((4*i).pairs, [(4,8)])
        self.assertEqual((4/i).pairs, [(2,4)])

        # Test for TypeError otherwise
        self.assertRaises(TypeError, lambda x : None + x, i)
        self.assertRaises(TypeError, lambda x : "foo" - x, i)
        self.assertRaises(TypeError, lambda x : set() * x, i)
        self.assertRaises(TypeError, lambda x : (lambda y : y) / x, i)


if __name__ == '__main__':
    unittest.main()
