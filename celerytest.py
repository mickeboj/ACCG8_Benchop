from cproj.tasks import solveproblem

PROBLEMS = ["prob1aI", "prob1bI", "prob1cI", "prob1aII","prob1bII", "prob1cII"]

if __name__ == '__main__':
    result = solveproblem.delay(PROBLEMS[0])
    res1 = result.get()
    print res1[0]    
