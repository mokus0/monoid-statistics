{-# LANGUAGE MultiParamTypeClasses #-}
module Math.Statistics.Monoid where

import qualified Math.Statistics as S
import Data.Monoid
import Debug.Trace

-- |these functions are more of an experiment than a usable library.
-- the stddev part in particular is cool but only numerically stable for values
-- with magnitude less than about 1e15 (for Double) and even then I haven't
-- tested it or analyzed it very thoroughly.  Mostly it's just a nifty
-- little generalization of one-pass stddev tricks to support arbitrary
-- concatenation of series summaries.

test = stddev (p1 `mappend` p2) - S.stddev (x++y)
    where
        w = 50000.0                                  
        p1 = mconcat (map pass1 x); x = [1..w]       
        p2 = mconcat (map pass1 y); y = [w+1..250000]

test2 = stddev (chunkReduce 1000 pass1 x) - S.stddev x
    where x = [1..10000]
test3 = stddev (chunkReduce 1000 pass1 $ map (1e10+) x) - S.stddev x
    where x = [1..10000]

-- test4 shows about where instability creeps into stddev algorithm for large values.
test4 = stddev (chunkReduce 10 pass1 $ map (+big) x) - S.stddev (map ((subtract big).(+big)) x)
    where 
        big = 1e17
        x = [1..10000]
test5 = (stddev (chunkReduce 1000 pass1 $ map (1e10+) x) - s) / s
    where
        x = (1e20:[1..10000])
        s = S.stddev x


chunkReduce :: Monoid b => Int -> (a -> b) -> [a] -> b
chunkReduce n f xs = chunkConcat n (map f xs)

chunkConcat :: Monoid a => Int -> [a] -> a
chunkConcat n xs = case map mconcat (chunk n xs) of
    [] -> mempty
    [x] -> x
    xs -> chunkConcat n xs

chunk :: Int -> [a] -> [[a]]
chunk n = go
    where 
        go []   = []
        go xs = c : go rest
            where (c,rest) = splitAt n xs


data Pass1 t = Pass1
    { p1sum     :: !t
    , p1count   :: !Int
    , p1a       :: !t   -- "excess" moment
    } deriving (Eq, Show)

instance RealFloat t => Monoid (Pass1 t) where
    mempty = Pass1 0 0 0
    mappend (Pass1 _ 0 _) p2 = p2
    mappend p1 (Pass1 _ 0 _) = p1
    mappend p1@(Pass1 s1 c1 a1) p2@(Pass1 s2 c2 a2) = Pass1 (s1 + s2) (c1+c2) (addM a1 a2)
        where
            addM = addExcessMoments s1 k1 s2 k2
            k1 = fromIntegral c1
            k2 = fromIntegral c2

expectedMoment mu q k = k*(mu-q)^2
addExcessMoments s1 k1 s2 k2 a1 a2 = a1 + a2 + u12
    where
        u12 = (k1*s2 - k2*s1)^2 / (k1*k2*(k1+k2))

-- only a useful monoid in the presence of a known mean:
data Pass2 t = Pass2
    { p2Moment2 :: !t
    , p2Moment3 :: !t
    , p2Moment4 :: !t
    } deriving (Eq, Show)

instance Num t => Monoid (Pass2 t) where
    mempty = Pass2 0 0 0
    mappend (Pass2 a1 b1 c1) (Pass2 a2 b2 c2) = Pass2 (a1+a2) (b1+b2) (c1+c2)

pass1 x = Pass1 x 1 0

pass2 p1 = p2
    where
        m = mean p1
        p2 x = Pass2 d2 d3 d4
            where
                d = x - m
                d2 = d*d
                d3 = d*d2
                d4 = d2*d2

-- pass1 stats
mean   (Pass1 s c _) = realToFrac s / fromIntegral c
var p1 = a / (k - 1)
    where
        a = p1a p1
        k = fromIntegral (p1count p1)
stddev p1 = sqrt (var p1)

x ~= y 
    = abs (x-y) < epsilon * max x y
    where epsilon = 1e-10


-- pass2 stats
moment2 (Pass1 _ c _) (Pass2 m2 _ _) = realToFrac m2 / fromIntegral c
moment3 (Pass1 _ c _) (Pass2 _ m3 _) = realToFrac m3 / fromIntegral c
moment4 (Pass1 _ c _) (Pass2 _ _ m4) = realToFrac m4 / fromIntegral c

skew p1 p2 = moment3 p1 p2 * moment2 p1 p2 ** (-1.5)
kurt p1 p2 = moment4 p1 p2 / (moment2 p1 p2 ^ 2) - 3

