from photometry import photometry
from simulation import simulation

class Image:
    def __init__(self, focus):
        self.focus = focus

        true_fwhm = 0.016*(focus-360)**2+6.0
        simulation(fwhm=true_fwhm)

        file_name = str(true_fwhm)
        file_name = file_name.replace('.', '_')
        self.fwhm = photometry(f'images/focus_{file_name}.fits')
        


def focus_finder(image1, image2, seen):

    print(f"image1: focus:{image1.focus}, FWHM:{image1.fwhm}")
    print(f"image2: focus:{image2.focus}, FWHM:{image2.fwhm} \n")

    if image1.fwhm > image2.fwhm:
        focus3 = image2.focus + (image2.focus - image1.focus)
        if focus3 in seen:
            return image2.focus
        seen.add(focus3)
        image3 = Image(focus3)
        return focus_finder(image2, image3, seen)
    elif image1.fwhm < image2.fwhm:
        focus3 = image1.focus - (image2.focus - image1.focus)
        if focus3 in seen:
            return image1.focus
        seen.add(focus3)
        image3 = Image(focus3)
        return focus_finder(image3, image1, seen)
    else:
        return image1.focus

    
import argparse

def main():
    parser = argparse.ArgumentParser(description="Find optimal focus for images based on FWHM.")
    parser.add_argument('--image1', type=int, default=380, help='Focus value for the first image')
    parser.add_argument('--image2', type=int, default=385, help='Focus value for the second image')
    args = parser.parse_args()
    # image1 = Image(335)
    # image2 = Image(340)
    # image1 = Image(380)
    # image2 = Image(385)

    image1 = Image(args.image1)
    image2 = Image(args.image2)

    seen = set()
    seen.add(image1.focus)
    seen.add(image2.focus)

    focus = focus_finder(image1, image2, seen)
    print(f"Optimal focus: {focus}")


if __name__ == "__main__":
    main()
