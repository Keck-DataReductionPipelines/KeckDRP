from KCWIPyDRP import kcwi_primitives
import argparse
import importlib



def go(iimg,rcp):
    if rcp is None:
        print("Checking %s header for recipe" % iimg)
        quit()

    mymodule = importlib.import_module("KCWIPyDRP.recipes."+str(rcp))
    myfunc = getattr(mymodule,rcp)

    p = kcwi_primitives.kcwi()

    print("executing recipe: %s" % rcp)
    myfunc(p)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=
        """Perform a reduction.
        """, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--recipe', type=str, help='reduction recipe')
    parser.add_argument('inimg', type=str, help='input image file')

    args = parser.parse_args()

    print("input  image: %s" % args.inimg)
    print("input recipe: %s" % args.recipe)

    go(args.inimg, args.recipe)

