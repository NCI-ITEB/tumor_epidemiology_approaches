/* Adapted from https://dakotaleemartinez.com/tutorials/how-to-add-active-highlight-to-table-of-contents/ */

document.addEventListener("DOMContentLoaded", (e) => {
      Scroller.init()
    });

function sleep(ms) {
    return new Promise(resolve => setTimeout(resolve, ms));
}

class Scroller {
  static init() {
    if(document.querySelector('.menu')) {
      this.contents = document.querySelectorAll('.menu a');
      this.contents.forEach(link => link.classList.add('transition', 'duration-200'))
      this.headers = Array.from(this.contents).map(link => {
        return document.querySelector(`#${link.href.split('#')[1]}`);
      })
      this.ticking = false;
      window.addEventListener('scroll', (e) => {
        this.onScroll()
        sleep(100)
      })
    }
  }

  static onScroll() {
    if(!this.ticking) {
      requestAnimationFrame(this.update.bind(this));
      this.ticking = true;
    }
  }

  static update() {
    this.activeHeader ||= this.headers[1];
    let activeIndex = this.headers.findIndex((header) => {
      return header.getBoundingClientRect().top > 180;
    });
    if(activeIndex == -1) {
      activeIndex = this.headers.length - 1;
    } else if(activeIndex > 0) {
      activeIndex--;
    }
    let active = this.headers[activeIndex];
    if(active !== this.activeHeader) {
      this.activeHeader = active;
      this.contents.forEach(link => link.classList.remove('is-active'));
      this.contents[activeIndex].classList.add('is-active');
    }
    this.ticking = false;
  }
}
