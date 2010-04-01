{-# LANGUAGE MultiParamTypeClasses #-}
-- |these functions are more of an experiment than a usable library.
-- the stddev part in particular is cool but only numerically stable for values
-- with magnitude less than about 1e15 (for Double) and even then I haven't
-- tested it or analyzed it very thoroughly.  Mostly it's just a nifty
-- little generalization of one-pass stddev tricks to support arbitrary
-- concatenation of series summaries.
module Math.Statistics.Monoid 
    ( Pass1(..), Pass2(..), Covar, AutoCovar
    , pass1, pass2
    , mean
    , var
    , stddev
    , moment2
    , moment3
    , moment4
    , skew
    , kurt
    , mkCovar, covar, correl, correlBy
    , covarA, covarB, covarAB
    , mkAutocovar, autoCovar, autoCorrel, autoCorrelBy
    , autoCovarPass1Stats
    ) where

-- import qualified Math.Statistics as S
import Data.Monoid

-- test = stddev (p1 `mappend` p2) - S.stddev (x++y)
--     where
--         w = 50000.0                                  
--         p1 = mconcat (map pass1 x); x = [1..w]       
--         p2 = mconcat (map pass1 y); y = [w+1..250000]
-- 
-- test2 = stddev (chunkReduce 1000 pass1 x) - S.stddev x
--     where x = [1..10000]
-- test3 = stddev (chunkReduce 1000 pass1 $ map (1e10+) x) - S.stddev x
--     where x = [1..10000]
-- 
-- -- test4 shows about where instability creeps into stddev algorithm for large values.
-- test4 = stddev (chunkReduce 10 pass1 $ map (+big) x) - S.stddev (map ((subtract big).(+big)) x)
--     where 
--         big = 1e17
--         x = [1..10000]
-- test5 = (stddev (chunkReduce 1000 pass1 $ map (1e10+) x) - s) / s
--     where
--         x = (1e20:[1..10000])
--         s = S.stddev x
-- 
-- 
-- chunkReduce :: Monoid b => Int -> (a -> b) -> [a] -> b
-- chunkReduce n f xs = chunkConcat n (map f xs)
-- 
-- chunkConcat :: Monoid a => Int -> [a] -> a
-- chunkConcat n xs = case map mconcat (chunk n xs) of
--     [] -> mempty
--     [x] -> x
--     xs -> chunkConcat n xs
-- 
-- chunk :: Int -> [a] -> [[a]]
-- chunk n = go
--     where 
--         go []   = []
--         go xs = c : go rest
--             where (c,rest) = splitAt n xs
-- 

data Pass1 t = Pass1
    { p1sum     :: !t
    , p1count   :: !Int
    , p1a       :: !t   -- "excess" moment
    , p1min     :: !t
    , p1max     :: !t
    } deriving (Eq, Show)

instance (Fractional t, Ord t) => Monoid (Pass1 t) where
    mempty = Pass1 0 0 0 0 0
    mappend (Pass1 _ 0 _ _ _) p2 = p2
    mappend p1 (Pass1 _ 0 _ _ _) = p1
    mappend p1@(Pass1 s1 c1 a1 mn1 mx1) p2@(Pass1 s2 c2 a2 mn2 mx2) = Pass1 (s1 + s2) (c1+c2) (addM a1 a2) (min mn1 mn2) (max mx1 mx2)
        where
            addM = addExcessMoments s1 k1 s2 k2
            k1 = fromIntegral c1
            k2 = fromIntegral c2

addExcessMoments s1 k1 s2 k2 a1 a2 = a1 + a2 + u12
    where
        u12 = (k1*s2 - k2*s1)^2 / (k1*k2*(k1+k2))

-- |For higher moments, 2 passes are required.  The 'Monoid' instance for 'Pass2'
-- is only valid when combining 'Pass2' objects that were constructed using
-- the same 'Pass1' object, and the final resulting values are only correct
-- when the final 'Pass2' object incorporates the same set of data as the
-- initial 'Pass1' object.
--
-- For example, to compute excess kurtosis:
--
-- > kurtFromScratch xs = kurt p1 (foldMap (pass2 p1) xs)
-- >    where p1 = foldMap pass1 xs
data Pass2 t = Pass2
    { p2Moment2 :: !t
    , p2Moment3 :: !t
    , p2Moment4 :: !t
    } deriving (Eq, Show)

instance Num t => Monoid (Pass2 t) where
    mempty = Pass2 0 0 0
    mappend (Pass2 a1 b1 c1) (Pass2 a2 b2 c2) = Pass2 (a1+a2) (b1+b2) (c1+c2)

pass1 x = Pass1 x 1 0 x x

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
mean   (Pass1 s c _ _ _) = s / fromIntegral c
var p1  | k > 1     = Just (a / (k - 1))
        | otherwise = Nothing
    where
        a = p1a p1
        k = fromIntegral (p1count p1)
stddev p1 = fmap (sqrt . realToFrac) (var p1)

-- pass2 stats
moment2 p1 (Pass2 m2 _ _) = m2 / fromIntegral (p1count p1)
moment3 p1 (Pass2 _ m3 _) = m3 / fromIntegral (p1count p1)
moment4 p1 (Pass2 _ _ m4) = m4 / fromIntegral (p1count p1)

skew p1 p2 = realToFrac (moment3 p1 p2) * realToFrac (moment2 p1 p2) ** (-1.5)
kurt p1 p2 = moment4 p1 p2 / (moment2 p1 p2 ^ 2) - 3

-- |1-pass covariance of 2 series.  Also incidentally computes all Pass1
-- statistics of both series as well as of the elementwise product of the 2 series.
data Covar t = Covar
    { covarA    :: !(Pass1 t)
    , covarB    :: !(Pass1 t)
    , covarAB   :: !(Pass1 t)
    } deriving (Eq, Show)

instance (Fractional t, Ord t) => Monoid (Covar t) where
    mempty = Covar mempty mempty mempty
    mappend (Covar a1 b1 ab1) (Covar a2 b2 ab2) 
        = Covar (mappend a1 a2) (mappend b1 b2) (mappend ab1 ab2)

-- |Convert a pair of corresponding elements from 2 different series into
-- the 'Covar' representation for the corresponding positions.
--
-- To build a 'Covar' for the whole series, do something like: 
-- > mconcat (zipWith mkCovar xs ys)
mkCovar a b = Covar (pass1 a) (pass1 b) (pass1 (a*b))

covar (Covar a b ab)
    | p1count a <= 1    = Nothing
    | otherwise         = Just $ (p1sum ab - mean b * p1sum a - mean a * p1sum b + p1sum a * mean b) / (fromIntegral (p1count a - 1))

correl c = correlBy id c

correlBy f c@(Covar a b ab) = do
    cov  <- covar c
    sd_a <- stddev a
    sd_b <- stddev b
    return (f cov / (sd_a * sd_b))

-- |1-pass autocovariance of a series with lag 1.
-- Also incidentally computes all other pass1 statistics of the series.
data AutoCovar t
    = EmptyAutoCovar
    | AutoCovar
        { autoCovarFirstA   :: !t
        , autoCovarLastB    :: !t
        , autoCovarMiddle   :: !(Covar t)
        } deriving (Eq, Show)

-- |Recover the Pass1 stats for the series from which
-- an 'AutoCovar' structure was built.
autoCovarPass1Stats (AutoCovar a _ c) = pass1 a `mappend` covarB c

mkAutocovar a = AutoCovar a a mempty

instance (Fractional a, Ord a) => Monoid (AutoCovar a) where
    mempty = EmptyAutoCovar
    mappend EmptyAutoCovar x = x
    mappend x EmptyAutoCovar = x
    mappend (AutoCovar firstA b0 c0) (AutoCovar b1 lastB c1)
        = AutoCovar firstA lastB (mconcat [c0, mkCovar b0 b1, c1])

autoCovar      (AutoCovar _ _ c) = covar c
autoCorrel     (AutoCovar _ _ c) = correl c
autoCorrelBy f (AutoCovar _ _ c) = correlBy f c