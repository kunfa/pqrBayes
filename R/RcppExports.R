# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

BBL <- function(xx, y, ystar, W, maxSteps, hatBeta, hatAlpha, hatInvTauSq, invSigAlpha0, hatLambdaSqStar, hatSigmaSq, aStar, bStar, alpha, gamma, progress) {
    .Call(`_pqrBayes_BBL`, xx, y, ystar, W, maxSteps, hatBeta, hatAlpha, hatInvTauSq, invSigAlpha0, hatLambdaSqStar, hatSigmaSq, aStar, bStar, alpha, gamma, progress)
}

BBLSS <- function(xx, y, ystar, W, maxSteps, hatAlpha, hatBeta, hatInvTauSq, invSigAlpha0, hatPi, hatLambdaSq, hatSigmaSq, aStar, bStar, alpha, gamma, sh1, sh0, progress) {
    .Call(`_pqrBayes_BBLSS`, xx, y, ystar, W, maxSteps, hatAlpha, hatBeta, hatInvTauSq, invSigAlpha0, hatPi, hatLambdaSq, hatSigmaSq, aStar, bStar, alpha, gamma, sh1, sh0, progress)
}

BGLPointMass <- function(xx, y, W, s, q, maxSteps, hatAlpha, hatBeta, hatInvTauSqStar, invSigAlpha0, hatPiStar, hatLambdaSqStar, hatSigmaSq, aStar, bStar, alpha, gamma, sh1, sh0, progress) {
    .Call(`_pqrBayes_BGLPointMass`, xx, y, W, s, q, maxSteps, hatAlpha, hatBeta, hatInvTauSqStar, invSigAlpha0, hatPiStar, hatLambdaSqStar, hatSigmaSq, aStar, bStar, alpha, gamma, sh1, sh0, progress)
}

BGL <- function(xx, y, W, s, q, maxSteps, hatBeta, hatAlpha, hatInvTauSq, invSigAlpha0, hatLambdaSqStar, hatSigmaSq, aStar, bStar, alpha, gamma, progress) {
    .Call(`_pqrBayes_BGL`, xx, y, W, s, q, maxSteps, hatBeta, hatAlpha, hatInvTauSq, invSigAlpha0, hatLambdaSqStar, hatSigmaSq, aStar, bStar, alpha, gamma, progress)
}

BL <- function(xx, y, W, maxSteps, hatBeta, hatAlpha, hatInvTauSq, invSigAlpha0, hatLambdaSqStar, hatSigmaSq, aStar, bStar, alpha, gamma, progress) {
    .Call(`_pqrBayes_BL`, xx, y, W, maxSteps, hatBeta, hatAlpha, hatInvTauSq, invSigAlpha0, hatLambdaSqStar, hatSigmaSq, aStar, bStar, alpha, gamma, progress)
}

BLSS <- function(xx, y, W, maxSteps, hatAlpha, hatBeta, hatInvTauSq, invSigAlpha0, hatPi, hatLambdaSq, hatSigmaSq, aStar, bStar, alpha, gamma, sh1, sh0, progress) {
    .Call(`_pqrBayes_BLSS`, xx, y, W, maxSteps, hatAlpha, hatBeta, hatInvTauSq, invSigAlpha0, hatPi, hatLambdaSq, hatSigmaSq, aStar, bStar, alpha, gamma, sh1, sh0, progress)
}

BRL <- function(xx, y, W, maxSteps, hatAlpha, hatBeta, hatTau, hatV, hatSg, invSigAlpha0, hatEtaSq, xi1, xi2, r1, a, b, progress) {
    .Call(`_pqrBayes_BRL`, xx, y, W, maxSteps, hatAlpha, hatBeta, hatTau, hatV, hatSg, invSigAlpha0, hatEtaSq, xi1, xi2, r1, a, b, progress)
}

BRLSS <- function(xx, y, W, maxSteps, hatAlpha, hatBeta, hatTau, hatV, hatSg, invSigAlpha0, hatPi, hatEtaSq, xi1, xi2, r1, a, b, sh1, sh0, progress) {
    .Call(`_pqrBayes_BRLSS`, xx, y, W, maxSteps, hatAlpha, hatBeta, hatTau, hatV, hatSg, invSigAlpha0, hatPi, hatEtaSq, xi1, xi2, r1, a, b, sh1, sh0, progress)
}

BRBL <- function(xx, y, ystar, W, maxSteps, hatAlpha, hatBeta, hatTau, hatV, hatSg, invSigAlpha0, hatEtaSq, xi1, xi2, r1, a, b, progress) {
    .Call(`_pqrBayes_BRBL`, xx, y, ystar, W, maxSteps, hatAlpha, hatBeta, hatTau, hatV, hatSg, invSigAlpha0, hatEtaSq, xi1, xi2, r1, a, b, progress)
}

BRBLSS <- function(xx, y, ystar, W, maxSteps, hatAlpha, hatBeta, hatTau, hatV, hatSg, invSigAlpha0, hatPi, hatEtaSq, xi1, xi2, r1, a, b, sh1, sh0, progress) {
    .Call(`_pqrBayes_BRBLSS`, xx, y, ystar, W, maxSteps, hatAlpha, hatBeta, hatTau, hatV, hatSg, invSigAlpha0, hatPi, hatEtaSq, xi1, xi2, r1, a, b, sh1, sh0, progress)
}

BRGL_SS <- function(xx, y, W, s, L, maxSteps, hatAlpha, hatBeta, hatTau, hatV, hatSg, invSigAlpha0, hatPi, hatEtaSq, xi1, xi2, r, a, b, sh1, sh0, progress) {
    .Call(`_pqrBayes_BRGL_SS`, xx, y, W, s, L, maxSteps, hatAlpha, hatBeta, hatTau, hatV, hatSg, invSigAlpha0, hatPi, hatEtaSq, xi1, xi2, r, a, b, sh1, sh0, progress)
}

BRGL <- function(xx, y, W, s, L, maxSteps, hatAlpha, hatBeta, hatTau, hatV, hatSg, invSigAlpha0, hatEtaSq, xi1, xi2, r, a, b, progress) {
    .Call(`_pqrBayes_BRGL`, xx, y, W, s, L, maxSteps, hatAlpha, hatBeta, hatTau, hatV, hatSg, invSigAlpha0, hatEtaSq, xi1, xi2, r, a, b, progress)
}

