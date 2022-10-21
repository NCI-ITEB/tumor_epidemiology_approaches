/*Adapted from https://www.aleksandrhovhannisyan.com/blog/how-to-add-a-copy-to-clipboard-button-to-your-jekyll-blog/*/
const copyCodeSnippetButtons = document.querySelectorAll('.copy-codeSnippet-button');
copyCodeSnippetButtons.forEach((copyCodeSnippetButton) => {
  const codeSnippet = copyCodeSnippetButton.previousSibling.innerText;

  copyCodeSnippetButton.addEventListener('click', () => {
    window.navigator.clipboard.writeText(codeSnippet);
    /*window.alert(code);*/
    copyCodeSnippetButton.classList.remove('fa-copy');
    copyCodeSnippetButton.classList.add('fa-check');

    setTimeout(() => {
      copyCodeSnippetButton.classList.remove('fa-check');
      copyCodeSnippetButton.classList.add('fa-copy');
    }, 2000);
  });
});
