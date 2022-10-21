/*Adapted from https://www.aleksandrhovhannisyan.com/blog/how-to-add-a-copy-to-clipboard-button-to-your-jekyll-blog/*/

const codeBlocks = document.querySelectorAll('.code-block-header + .highlighter-rouge');
const copyCodeBlockButtons = document.querySelectorAll('.copy-codeBlock-button');

copyCodeBlockButtons.forEach((copyCodeBlockButton, index) => {
  const code = codeBlocks[index].innerText;

  copyCodeBlockButton.addEventListener('click', () => {
    window.navigator.clipboard.writeText(code);
    copyCodeBlockButton.classList.remove('fa-copy');
    copyCodeBlockButton.classList.add('fa-check');

    setTimeout(() => {
      copyCodeBlockButton.classList.remove('fa-check');
      copyCodeBlockButton.classList.add('fa-copy');
    }, 2000);
  });
});
