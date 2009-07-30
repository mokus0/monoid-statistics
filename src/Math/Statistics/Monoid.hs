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
    , p1a       :: !t   -- "excess" moment about 0
    , p1b       :: !t   -- "excess" moment about Q
    , p1q       :: t    -- Q
    } deriving (Eq, Show)

instance RealFloat t => Monoid (Pass1 t) where
    mempty = Pass1 0 0 0 0 1
    mappend (Pass1 _ 0 _ _ _) p2 = p2
    mappend p1 (Pass1 _ 0 _ _ _) = p1
    mappend p1@(Pass1 s1 c1 a1 b1 q1) p2@(Pass1 s2 c2 a2 b2 q2) = Pass1 (s1 + s2) (c1+c2) (addM a1 a2) (addM b1' b2') $! q
        where
            q = balanceQ p1 p2
            b1' = recenterQ q p1
            b2' = recenterQ q p2
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

-- |"Magic" constant used to compute variance in pass 1 by accumulation of
-- the values p1a and p1b (moments about 0 and q respectively).  The further q
-- is from the actual mean of the set, the more truncation error will be
-- introduced in the computation for variance and stddev.
-- In practice, the stability of this algorithm seems to be almost totally 
-- insensitive to the value returned, so currently it just returns a
-- constant to avoid wasted time and loss of precision from too many 
-- recenters.
q :: (Fractional a, Ord a) => a -> a
q x = 5000.5 -- head (dropWhile ((< x)) qs)
    where
        qs = iterate (*16) 16

pass1 x = Pass1 x 1 0 0 qx
    where
        qx = q x

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
mean   (Pass1 s c _    _    _) = realToFrac s / fromIntegral c
var p1 = recenterQ (mean p1) p1 / (k - 1)
    where
        k = fromIntegral (p1count p1)
stddev p1 = sqrt (var p1)

center p1 = recenter (mean p1) p1

recenter  newQ p1@(Pass1 s c a_2r b_2r q)
    = Pass1 s c a_2r (recenterQ newQ p1) newQ

recenterQ newQ p1@(Pass1 s c a_2r b_2r q)
    | q == newQ
    = b_2r
recenterQ newQ p1@(Pass1 _ 0 _ _ _) = 0
recenterQ newQ p1@(Pass1 _ 1 _ _ _) = 0
recenterQ newQ p1@(Pass1 _ _ a b q)
    | b == a    = a
    | otherwise
    = (b - a) * (newQ / q) + a
            

x ~= y 
    = abs (x-y) < epsilon * max x y
    where epsilon = 1e-10


balanceQ (Pass1 _ 1 _ _ _) (Pass1 _ _ _ _ q) = q
balanceQ (Pass1 _ _ _ _ q) (Pass1 _ 1 _ _ _) = q
balanceQ p1@(Pass1 s1 c1 _ _ q1) p2@(Pass1 s2 c2 _ _ q2)
    -- | q1 == q2 |
    | q1 ~= q2
    = max q1 q2
    
    | otherwise
    = q $ (mean2 p1 p2)
    
    where
        mean2 (Pass1 s1 c1 _ _ _) (Pass1 s2 c2 _ _ _) = (s1 + s2) / fromIntegral (c1+c2)
        x ~= y
            = max x y <= alpha * min x y
            where alpha = 2
            -- = abs (x-y) < epsilon * max x y
            -- where epsilon = 1e1

-- pass2 stats
moment2 (Pass1 _ c _ _ _) (Pass2 m2 _ _) = realToFrac m2 / fromIntegral c
moment3 (Pass1 _ c _ _ _) (Pass2 _ m3 _) = realToFrac m3 / fromIntegral c
moment4 (Pass1 _ c _ _ _) (Pass2 _ _ m4) = realToFrac m4 / fromIntegral c

skew p1 p2 = moment3 p1 p2 * moment2 p1 p2 ** (-1.5)
kurt p1 p2 = moment4 p1 p2 / (moment2 p1 p2 ^ 2) - 3

